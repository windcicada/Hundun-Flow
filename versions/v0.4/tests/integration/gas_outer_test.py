#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Fixed reference sweeps preserve physical audits and restart identity."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary, mpi = [str(Path(p).resolve()) for p in sys.argv[1:3]]
source = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location('evidence', source/'tools/v04_evidence_validate.py')
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)

def call(n, args, okay=True):
    result = subprocess.run([mpi, '--oversubscribe', '--bind-to', 'none', '-n', str(n),
                             binary, *map(str, args)], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, universal_newlines=True, timeout=60)
    assert (result.returncode == 0) == okay, result.stdout[-10000:]
    return result.stdout

with tempfile.TemporaryDirectory(prefix='hf-outer-') as tmp:
    root = Path(tmp)
    case = root/'c'
    call(1, ['init-case', '--output', case])
    model = json.loads((case/'case.json').read_text())
    model['mesh']['exact_cells'] = [4, 4, 4]
    model['mesh']['minimum_spacing'] = [.25]*3
    model['flow']['pressure_reference'] = 'boundary_absolute'
    for face in model['boundaries'].values():
        face['flow_kind'] = 'periodic'
        face['thermal_kind'] = 'none'
        for key, value in list(face.items()):
            if isinstance(value, float): face[key] = 0.0
            elif isinstance(value, list) and len(value) == 3: face[key] = [0., 0., 0.]
        face['mach_limit'] = .95
    model['boundaries']['x_min'].update(flow_kind='velocity_inlet', velocity=[1., 0., 0.], temperature=300.)
    model['boundaries']['x_max'].update(flow_kind='pressure_outlet', pressure=101325., backflow_temperature=300.)
    model['time'].update(control='fixed', scheme='cn_be', initial_dt=.001, maximum_dt=.001,
                         convective_cfl_definition='directional_max')
    model['solver']['coupling'] = 'outer_corrected'
    def write(): (case/'case.json').write_text(json.dumps(model))
    def run(name, n=2, restart=None, okay=True):
        write()
        args = ['run', case, '--output', root/name, '--steps', 1,
                '--output-interval', 0, '--restart-interval', 1]
        args += ['--restart', restart] if restart else ['--initial-state', '101325,300,1,0,0']
        text = call(n, args, okay)
        if not okay: return text
        row = json.loads((root/name/'evidence.jsonl').read_text().splitlines()[-1])
        validator.validate_v6_v8_runtime_record(row, 1, 8)
        expected = model['solver'].get('reference_outer_iterations', 0)
        assert row['cold']['reference_outer_iterations'] == expected
        if expected:
            assert row['cold']['outer_iterations'] == expected
            assert text.count('cold_outer outer=') == expected
        return row
    natural = run('natural')
    for count in (1, 3):
        model['solver']['reference_outer_iterations'] = count
        write()
        assert 'reference_outer_iterations='+str(count) in call(2, ['check', case])
        row = run('n'+str(count))
        bad = copy.deepcopy(row)
        bad['cold']['reference_outer_iterations'] = count+1
        try: validator.validate_v6_v8_runtime_record(bad, 1, 8)
        except validator.EvidenceError: pass
        else: raise AssertionError('forged fixed schedule admitted')
    resumed = run('resume', n=4, restart=root/'n3/Restart')
    assert resumed['step'] == 2
    model['solver']['reference_outer_iterations'] = 1
    run('drift', restart=root/'n3/Restart', okay=False)
    # A genuine thermal transient cannot meet near-roundoff original-equation
    # tolerances in one sweep. Fixed scheduling must preserve the rejected step.
    model['boundaries']['x_min']['temperature'] = 310.
    model['time']['minimum_dt'] = .001
    model['solver']['cold_stopping'] = dict(reference_time=1., momentum=1e-25,
                                          enthalpy=1e-25, species=1e-25)
    failure = run('residual', okay=False)
    assert 'detail=17837' in failure, failure[-8000:]
    assert failure.count('cold_outer outer=') == 1, failure[-8000:]
    assert 'cold_terminal_audit' in failure and 'passed=0' in failure, failure[-8000:]
    assert not (root/'residual/Restart/current').exists()
    evidence = root/'residual/evidence.jsonl'
    assert not evidence.exists() or not evidence.read_text().strip()
    model['solver'].pop('cold_stopping')
    for invalid in (-1, 65, 1.5, True):
        model['solver']['reference_outer_iterations'] = invalid
        write()
        call(2, ['check', case], okay=False)
    model['solver']['reference_outer_iterations'] = 1
    model['time']['scheme'] = 'backward_euler'
    model['solver']['coupling'] = 'PISO'
    write()
    call(2, ['check', case], okay=False)
    print('Fixed outer: 1/3 sweeps, residual audit rejection, MPI 2-to-4 restart, identity drift, strict input and evidence PASS')
