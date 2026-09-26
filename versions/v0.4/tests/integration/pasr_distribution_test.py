#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Spatial PaSR sources and native recovery agree on 1/2/4 MPI ranks."""
import io
import json
import math
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import tempfile

binary, launcher = sys.argv[1:3]
rank_flag = sys.argv[3] if len(sys.argv) > 3 else '-n'
pre_flags = [v for v in sys.argv[4].split(';') if v] if len(sys.argv) > 4 else []
post_flags = [v for v in sys.argv[5].split(';') if v] if len(sys.argv) > 5 else []


def call(ranks, args):
    command = [launcher, rank_flag, str(ranks)] + pre_flags + [binary] + post_flags + list(map(str, args))
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True, timeout=60)
    assert result.returncode == 0, result.stdout[-10000:]


def frame(root, step):
    owned = {}
    for path in (root/'Visit').glob('s%dr*.vtk' % step):
        stream = io.BytesIO(path.read_bytes())
        def line():
            return stream.readline().decode().strip()
        assert line() == '# vtk DataFile Version 3.0'
        line()
        assert line() == 'BINARY' and line() == 'DATASET STRUCTURED_GRID'
        assert line().startswith('DIMENSIONS ')
        words = line().split()
        assert words[0] == 'POINTS' and words[2] == 'double'
        count = int(words[1])
        xyz = list(struct.iter_unpack('>3d', stream.read(count*24)))
        assert line() == '' and line() == 'POINT_DATA '+str(count)
        fields = {}
        while stream.tell() < len(stream.getbuffer()):
            words = line().split()
            if not words:
                continue
            if words[0] == 'FIELD':
                name, width, points, dtype = line().split()
                assert int(points) == count
            elif words[0] == 'VECTORS':
                _, name, dtype = words
                width = 3
            else:
                _, name, dtype, width = words
                assert line() == 'LOOKUP_TABLE default'
            width = int(width)
            code = 'd' if dtype == 'double' else 'B'
            fields[name] = list(struct.iter_unpack('>'+code*width,
                stream.read(count*width*struct.calcsize(code))))
            assert line() == ''
        for i, point in enumerate(xyz):
            if fields['vtkGhostType'][i][0]:
                continue
            assert point not in owned
            owned[point] = {name: rows[i] for name, rows in fields.items()
                            if name != 'vtkGhostType'}
    assert len(owned) == 512
    return owned


def compare(left, right):
    assert left.keys() == right.keys()
    for cell in left:
        assert left[cell].keys() == right[cell].keys()
        for field in left[cell]:
            for a, b in zip(left[cell][field], right[cell][field]):
                assert math.isfinite(a) and math.isfinite(b)
                assert abs(a-b) <= 1e-9*max(1., abs(a), abs(b)), (cell, field, a, b)


with tempfile.TemporaryDirectory(prefix='hundun-pasr-distribution-') as directory:
    work = Path(directory)
    case = work/'case'
    shutil.copytree(str(Path(__file__).resolve().parents[1]/'fixtures'/'reacting-cantera-isomer'), str(case))
    model = json.loads((case/'case.json').read_text())
    model['mesh']['domain']['upper'] = [.008]*3
    model['mesh']['minimum_spacing'] = [.001]*3
    model['time'].update(scheme='cn_be', initial_dt=1e-6, minimum_dt=1e-6, maximum_dt=1e-6)
    model['solver']['coupling'] = 'outer_corrected'
    model['solver']['cold_stopping'] = dict(reference_time=1e-6, momentum=1e-8,
                                            enthalpy=1e-9, species=1e-9)
    model['flow']['pressure_reference'] = 'boundary_absolute'
    for face in ('z_min', 'z_max'):
        model['boundaries'][face].update(flow_kind='pressure_outlet', pressure=101325.,
            temperature=500., total_temperature=500., backflow_temperature=500.,
            allow_backflow=True, scalars=[dict(stable_name='A', kind='zero_gradient',
                value=.5, backflow_kind='dirichlet', backflow_value=.5)])
    model['reaction'].update(model='auto', ensemble={'fields': 1},
                             mixing={'c_z': 1., 'turbulent_schmidt': .7})
    (case/'case.json').write_text(json.dumps(model))
    thermo = (case/'thermophysics.d').read_text()
    properties = []
    for name in ('A', 'B'):
        block = thermo.split('species '+name+'\n')[1].split('end_species')[0]
        mw = float(re.search(r'molecular_weight (\S+)', block).group(1))
        nasa = list(map(float, re.search(r'nasa7_low (.*)', block).group(1).split()))
        properties.append((8314.46261815324/mw, nasa))
    transfer = work/'transfer'
    transfer.mkdir()
    (transfer/'state.txt').write_text('HUNDUN_PDF_TRANSFER 1\n8 8 8 17 .001 1e-6 101325 1 2\nA B\n')
    flow, pdf, rho = [], [], []
    for z in range(8):
        for y in range(8):
            for x in range(8):
                a = .5 + .4*math.sin(2*math.pi*(x+.5)/8)*math.cos(2*math.pi*(z+.5)/8)
                temperature = 500 + 100*math.sin(2*math.pi*(y+.5)/8)
                h = gas = 0.
                for fraction, (r, c) in zip((a, 1-a), properties):
                    gas += fraction*r
                    h += fraction*r*(sum(c[i]*temperature**(i+1)/(i+1) for i in range(5))+c[5])
                flow.extend((1., 0., 0., 101325.))
                pdf.extend((a, 1-a, h))
                rho.append(101325/(gas*temperature))
    for name, values in [('flow.f64', flow), ('pdf0.f64', pdf), ('rho_ref.f64', rho)]:
        (transfer/name).write_bytes(struct.pack('<'+'d'*len(values), *values))
    (transfer/'fluid.u8').write_bytes(bytes([1])*512)
    seed = work/'seed'
    call(1, ['import', transfer, '--format', 'pdf-transfer-v1', '--case', case, '--output', seed])
    runs = {}
    for ranks in (1, 2, 4):
        output = work/('r%d' % ranks)
        call(ranks, ['run', case, '--restart', seed, '--output', output, '--steps', '2',
                      '--output-interval', '1', '--restart-interval', '1', '--diagnostics-interval', '1'])
        runs[ranks] = frame(output, 19)
        initial_a = sum(r*pdf[3*i] for i, r in enumerate(rho))*1e-9
        final_a = sum(v['Density'][0]*v['A'][0] for v in runs[ranks].values())*1e-9
        assert final_a < initial_a*(1-1e-8), (initial_a, final_a)
        for row in map(json.loads, (output/'diagnostics.jsonl').read_text().splitlines()):
            budget = row['payload']
            assert abs(budget['mass_balance_defect_kg_s']*budget['dt']) < 1e-10*budget['mass_kg']
            scale = max(1., abs(budget['internal_energy_J']) + budget['kinetic_energy_J'])
            assert abs(budget['total_energy_balance_defect_W']*budget['dt']) < 1e-10*scale
            composition = budget['composition_balance']
            assert any(abs(s['chemistry_source']) > 0 for s in composition['species'])
            for entry in composition['species'] + composition['elements']:
                assert entry['relative_defect'] < 1e-6, entry
    compare(runs[1], runs[2])
    compare(runs[1], runs[4])
    reference, resumed = work/'reference', work/'resumed'
    for ranks, source, output in [(1, work/'r1', reference), (4, work/'r2', resumed)]:
        call(ranks, ['run', case, '--restart', source/'Restart', '--output', output,
                      '--steps', '1', '--output-interval', '1', '--restart-interval', '0'])
    compare(frame(reference, 20), frame(resumed, 20))
print('PASS: spatial real-Cantera PaSR fields, budgets and 2-to-4 recovery')
