#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""CN/BE passive transport: conservative rows, affine invariance and Restart."""
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

with tempfile.TemporaryDirectory(prefix='hf-passive-') as tmp:
    root = Path(tmp)
    case = root/'c'
    call(1, ['init-case', '--output', case])
    model = json.loads((case/'case.json').read_text())
    model['mesh']['exact_cells'] = [4, 4, 4]
    model['mesh']['minimum_spacing'] = [.25]*3
    model['flow'].update(pressure_reference='boundary_absolute')
    model['transported_scalars'] = [dict(stable_name=name, role='passive_scalar',
        molecular_schmidt=1., turbulent_schmidt=1.) for name in ('mix', 'signed')]
    for face in model['boundaries'].values():
        face['flow_kind'] = 'periodic'
        face['thermal_kind'] = 'none'
        for key, value in list(face.items()):
            if isinstance(value, float): face[key] = 0.0
            elif isinstance(value, list) and len(value) == 3: face[key] = [0., 0., 0.]
        face.update(mach_limit=.95, relaxation=1., scalars=[])
    for name, inlet in [('x_min', True), ('x_max', False)]:
        face=model['boundaries'][name]
        face.update(flow_kind='velocity_inlet' if inlet else 'pressure_outlet',
            velocity=[1.,0.,0.] if inlet else [0.,0.,0.], temperature=310.,
            pressure=101325., backflow_temperature=310.,
            scalars=[dict(stable_name=name, kind='dirichlet' if inlet else 'zero_gradient',
                value=value if inlet else 0., backflow_kind='dirichlet', backflow_value=value)
                for name,value in [('mix',.8),('signed',.4)]])
    dt=.001
    model['time'].update(control='fixed', scheme='cn_be', initial_dt=dt,
        maximum_dt=dt, minimum_dt=dt)
    model['solver']['coupling']='outer_corrected'
    (case/'case.json').write_text(json.dumps(model))
    assert 'passive_workspace_bytes=0 ' not in call(2,['check',case])
    def run(name,n,steps,restart=None,initial='101325,300,1,0,0,0.25,-1.25'):
        args=['run',case,'--output',root/name,'--steps',steps,
            '--output-interval',1,'--restart-interval',1]
        args+=['--restart',restart] if restart else ['--initial-state',initial]
        text=call(n,args)
        rows=[json.loads(line) for line in (root/name/'evidence.jsonl').read_text().splitlines()]
        monitors={json.loads(line)['step']:json.loads(line)['payload']
            for line in (root/name/'monitor.jsonl').read_text().splitlines()}
        for row in rows:
            validator.validate_v6_v8_runtime_record(row,1,8)
            cold=row['cold']
            assert cold['passive_scalar_count']==2 and cold['passive_iterations']>0
            assert cold['passive_residual']<128*sys.float_info.epsilon
            assert cold['passive_balance_defect']<128*sys.float_info.epsilon
            for key in ('passive_scalar_count','passive_solve_calls','passive_iterations',
                        'passive_residual','passive_balance_defect'):
                assert monitors[row['step']][key]==cold[key],key
        return rows[-1]
    last=run('all',2,2)
    run('first',2,1)
    resumed=run('next',4,1,root/'first/Restart')
    assert resumed['step']==last['step']==2
    def samples(name):
        vtk=sorted((root/name/'Visit').glob('*.vtk'))
        # The index selects the last frame, retaining every block.
        lines=(root/name/'Visit/solution.visit').read_text().splitlines()
        n=int(lines[0].split()[1]); names=lines[-n:]
        result={}
        for filename in names:
            data=fields(root/name/'Visit'/filename)
            for i in range(len(data['Points'])//3):
                xyz=tuple(data['Points'][3*i:3*i+3])
                if not all(0.<x<1. for x in xyz): continue
                values={key:tuple(value[3*i:3*i+3]) if key=='Velocity' else value[i]
                    for key,value in data.items() if key!='Points'}
                result[xyz]=values
        return result
    a,b=samples('all'),samples('next')
    assert len(a)==len(b)==64
    max_affine=max_restart=0.
    for xyz,left in a.items():
        right=b[xyz]
        q=left['mix']; s=left['signed']
        assert .25-1e-12<=q<=.8+1e-12, (xyz,q)
        max_affine=max(max_affine,abs(s-(3*q-2)))
        for key in ('mix','signed','Density','Temperature','Enthalpy','PressureGauge'):
            error=abs(left[key]-right[key])/max(1.,abs(left[key]),abs(right[key]))
            max_restart=max(max_restart,error)
    assert max_affine<1e-10 and max_restart<1e-9, (max_affine,max_restart)
    assert max(row['mix'] for row in a.values())>.251
    assert max(abs(row['Density']-1.1766345660921333) for row in a.values())>1e-6
    forged=copy.deepcopy(resumed); forged['cold']['passive_residual']=1e-4
    try: validator.validate_v6_v8_runtime_record(forged,1,8)
    except validator.EvidenceError: pass
    else: raise AssertionError('forged passive residual admitted')
    print('CN/BE passive: variable density, signed affine tracer, 2-to-4 restart PASS',
        dict(affine=max_affine,restart=max_restart,residual=last['cold']['passive_residual'],
             balance=last['cold']['passive_balance_defect']))

    # An internal cube keeps solid scalar states while fluid-side diffusion
    # transports both tracers through MPI interfaces.
    model['mesh'].update(domain=dict(lower=[-2.]*3, upper=[2.]*3), exact_cells=[12]*3,
        minimum_spacing=[1./3]*3, immersed_boundary=dict(stl_file='cube.stl',fluid_side='outside'))
    model['mesh']['limits']['max_global_cells']=1728
    model['mesh']['limits']['max_memory_bytes_per_rank']=1073741824
    (case/'cube.stl').write_bytes((source/'versions/v0.4/tests/data/cube_binary.stl').read_bytes())
    for boundary in model['boundaries'].values():
        if boundary['flow_kind']!='periodic':
            boundary.update(temperature=300.,backflow_temperature=300.,
                velocity=[1e-6,0.,0.] if boundary['flow_kind']=='velocity_inlet' else [0.,0.,0.])
    (case/'case.json').write_text(json.dumps(model))
    ibm=run('ibm',2,1,initial='101325,300,0,0,0,0.25,-1.25')
    points={}
    for path in (root/'ibm/Visit').glob('*.vtk'):
        data=fields(path)
        for i in range(len(data['Points'])//3):
            xyz=tuple(data['Points'][3*i:3*i+3])
            if all(-2.<x<2. for x in xyz): points[xyz]=(data['mix'][i],data['signed'][i])
    solids=[v for xyz,v in points.items() if all(abs(x)<1. for x in xyz)]
    assert len(solids)==216 and all(v==(.25,-1.25) for v in solids),len(solids)
    assert max(v[0] for v in points.values())>.25+1e-9
    assert max(abs(v[1]-(3*v[0]-2)) for v in points.values())<1e-10
    print('IBM passive: 216 solid values held, fluid diffusion and affine closure PASS',
        ibm['cold']['passive_residual'])
