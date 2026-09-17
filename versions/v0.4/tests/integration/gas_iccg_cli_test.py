#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Native IC/PCG: original residuals, SPD admission and cross-rank Restart."""
import copy
import importlib.util
import json
import re
import struct
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
    result = subprocess.run([mpi, '--oversubscribe', '--bind-to', 'none', '-n', str(n), binary,
                             *map(str, args)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True, timeout=90)
    assert (result.returncode == 0) == okay, result.stdout[-12000:]
    return result.stdout

def fields(path):
    with path.open('rb') as stream:
        for _ in range(5): stream.readline()
        point_header = stream.readline().split()
        assert point_header[0] == b'POINTS' and point_header[2] == b'double', point_header
        count = int(point_header[1])
        points = struct.unpack('>'+str(3*count)+'d', stream.read(24*count))
        assert stream.readline() == b'\n'
        assert stream.readline().split() == [b'POINT_DATA', str(count).encode()]
        result = {'Points': points}
        while True:
            header = stream.readline().split()
            if not header: break
            name = header[1].decode()
            width = 3 if header[0] == b'VECTORS' else int(header[3])
            if header[0] == b'SCALARS': assert stream.readline().split() == [b'LOOKUP_TABLE', b'default']
            code, size = ('d', 8) if header[2] == b'double' else ('B', 1)
            result[name] = struct.unpack('>'+str(width*count)+code, stream.read(size*width*count))
            assert stream.readline() == b'\n'
        return result

with tempfile.TemporaryDirectory(prefix='hf-iccg-') as tmp:
    root = Path(tmp)
    case = root/'c'
    call(1, ['init-case', '--output', case])
    model = json.loads((case/'case.json').read_text())
    model['mesh']['exact_cells'] = [4, 4, 4]
    model['mesh']['minimum_spacing'] = [.25]*3
    model['flow'].update(pressure_reference='boundary_absolute', thermodynamic_pressure_pa=101325.)
    for face in model['boundaries'].values():
        face['flow_kind'] = 'periodic'
        face['thermal_kind'] = 'none'
        for key, value in list(face.items()):
            if isinstance(value, float): face[key] = 0.0
            elif isinstance(value, list) and len(value) == 3: face[key] = [0., 0., 0.]
        face['mach_limit'] = .95
    model['boundaries']['x_min'].update(flow_kind='velocity_inlet', velocity=[.125, 0., 0.], temperature=300.)
    model['boundaries']['x_max'].update(flow_kind='pressure_outlet', pressure=101325., backflow_temperature=300.)
    dt = 1./1024.
    model['time'].update(control='fixed', scheme='cn_be', initial_dt=dt, maximum_dt=dt,
                         minimum_dt=dt, convective_cfl_definition='directional_max')
    model['solver']['coupling'] = 'outer_corrected'
    model['solver']['pressure_linear'].update(algorithm='pcg', krylov_restart=0,
                                       absolute_tolerance=1e-13, relative_tolerance=1e-13)
    logs = {}
    def write(): (case/'case.json').write_text(json.dumps(model))
    def run(name, n=2, restart=None, okay=True, initial='101325,300,0,0,0'):
        write()
        args = ['run', case, '--output', root/name, '--steps', 1,
                '--output-interval', 1, '--restart-interval', 1, '--observe-mg-cost']
        args += ['--restart', restart] if restart else ['--initial-state', initial]
        text = call(n, args, okay)
        logs[name] = text
        if not okay: return text
        row = json.loads((root/name/'evidence.jsonl').read_text().splitlines()[-1])
        validator.validate_v6_v8_runtime_record(row, 1, 8)
        cold = row['cold']
        iccg = model['solver']['pressure_linear']['algorithm'] == 'pcg'
        assert cold['pressure_iccg'] == iccg
        assert cold['pressure_original_l2'] <= cold['pressure_original_l2_limit']
        if iccg: assert 'kind=iccg' in text and 'spd_admitted=1' in text
        monitor = json.loads((root/name/'monitor.jsonl').read_text().splitlines()[-1])
        # Monitor wraps its payload in the stable output envelope.
        assert 'pressure_iccg' in json.dumps(monitor)
        return row
    write()
    assert 'pressure_algorithm=pcg_iccg_spd' in call(2, ['check', case])
    fresh = run('new')
    assert fresh['cold']['pressure_iterations'] > 0
    resumed = run('next', n=4, restart=root/'new/Restart')
    assert resumed['step'] == 2
    forged = copy.deepcopy(resumed)
    forged['cold']['pressure_original_l2'] = 2*forged['cold']['pressure_original_l2_limit']
    try: validator.validate_v6_v8_runtime_record(forged, 1, 8)
    except validator.EvidenceError: pass
    else: raise AssertionError('ICCG original residual gate forged')
    model['solver']['pressure_linear'].update(algorithm='fgmres', krylov_restart=12)
    baseline = run('fgmres')
    matrix = r'cold_matrix outer=0 matrix=([0-9a-f]+) rhs=([0-9a-f]+)'
    assert re.search(matrix, logs['new']).groups() == re.search(matrix, logs['fgmres']).groups()
    differences = {}
    for path in sorted((root/'new/Visit').glob('*.vtk')):
        left, right = fields(path), fields(root/'fgmres/Visit'/path.name)
        for name in ('Points', 'Velocity', 'Temperature', 'PressureGauge', 'Density', 'Enthalpy'):
            a, b = left[name], right[name]
            assert len(a) == len(b)
            error = max(abs(x-y) for x, y in zip(a, b))/max(1., max(map(abs, a)), max(map(abs, b)))
            differences[name] = max(differences.get(name, 0.), error)
            assert error <= 1e-11, (name, error)
    assert differences, 'paired field output missing'
    print('ICCG/FGMRES same first matrix/RHS and step fields:', differences)
    model['solver']['pressure_linear'].update(algorithm='pcg', krylov_restart=0)
    model['flow'].pop('thermodynamic_pressure_pa')
    rejected = run('nonsym', okay=False, initial='101325,300,.125,0,0')
    assert '17871' in rejected, rejected[-10000:]
    assert not (root/'nonsym/Restart/current').exists()
    model['time']['scheme'] = 'backward_euler'
    model['solver']['coupling'] = 'PISO'
    write()
    call(2, ['check', case], okay=False)
    print('Native ICCG: SPD pressure, original L2/continuity, nonsymmetric rejection and MPI 2-to-4 restart PASS')
