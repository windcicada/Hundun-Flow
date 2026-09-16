#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Distinguish directional/outgoing admission using a diagonal uniform flow."""
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
validator = importlib.util.module_from_spec(spec); spec.loader.exec_module(validator)
def call(n,args,okay=True):
    p=subprocess.run([mpi,'--oversubscribe','--bind-to','none','-n',str(n),binary,*map(str,args)],
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True,timeout=60)
    assert (p.returncode==0)==okay,p.stdout[-8000:]
    return p.stdout
with tempfile.TemporaryDirectory(prefix='hf-cfl-') as tmp:
    root=Path(tmp);case=root/'c'
    call(1,['init-case','--output',case])
    model=json.loads((case/'case.json').read_text())
    model['mesh']['exact_cells']=[4,4,4];model['mesh']['minimum_spacing']=[.25]*3
    model['flow']['pressure_reference']='boundary_absolute'
    for face in model['boundaries'].values():
        face['flow_kind']='periodic';face['thermal_kind']='none'
        for key,value in list(face.items()):
            if isinstance(value,float):face[key]=0.0
            elif isinstance(value,list) and len(value)==3:face[key]=[0.,0.,0.]
        face['mach_limit']=.95
    model['boundaries']['x_min'].update(flow_kind='velocity_inlet',velocity=[1.,1.,1.],temperature=300.)
    model['boundaries']['x_max'].update(flow_kind='pressure_outlet',pressure=101325.,backflow_temperature=300.)
    model['time'].update(control='fixed',initial_dt=.075,maximum_dt=.075,
        convective_cfl=.3,convective_cfl_margin=.05,convective_cfl_definition='directional_max')
    def run(name,steps=1,restart=None,okay=True,n=2):
        (case/'case.json').write_text(json.dumps(model))
        if okay:
            assert 'cfl_definition='+model['time']['convective_cfl_definition'] in call(n,['check',case])
        args=['run',case,'--output',root/name,'--steps',steps,'--output-interval',0,'--restart-interval',1]
        args+=['--restart',restart] if restart else ['--initial-state','101325,300,1,1,1']
        text=call(n,args,okay)
        if not okay:return text
        rows=[json.loads(s) for s in (root/name/'evidence.jsonl').read_text().splitlines()]
        for row in rows:validator.validate_v6_v8_runtime_record(row,1,8)
        return rows
    for scheme in ('cn_be','backward_euler'):
        model['time']['scheme']=scheme
        model['solver']['coupling']='PISO'
        model['time']['convective_cfl_definition']='directional_max'
        rows=run(scheme)
        audit=rows[-1]['terminal_physical_audit']['committed_convective_cfl']
        assert abs(audit['directional_max']-.3)<1e-11,audit
        assert abs(audit['out_max']-.9)<1e-11,audit
        assert audit['definition']=='directional_max'
        corrupted=copy.deepcopy(rows[-1]);corrupted['terminal_physical_audit']['committed_convective_cfl']['directional_winner']['maximum_face_mass_flow']*=2
        try:validator.validate_v6_v8_runtime_record(corrupted,1,8)
        except validator.EvidenceError:pass
        else:raise AssertionError('forged directional winner admitted')
        resumed=run(scheme+'r',restart=root/scheme/'Restart',n=4)
        assert resumed[-1]['step']==2
        restored=resumed[-1]['terminal_physical_audit']['committed_convective_cfl']
        assert restored['definition']=='directional_max' and abs(restored['directional_max']-.3)<1e-11
        model['time']['control']='adaptive_flow'
        direct=run(scheme+'direct')[-1]['terminal_physical_audit']['committed_convective_cfl']
        assert abs(direct['dt']-.075)<1e-12,direct
        model['time']['control']='fixed'
        model['time']['convective_cfl_definition']='outgoing_sum'
        assert 'detail=10211' in run(scheme+'bad',okay=False)
        run(scheme+'drift',restart=root/scheme/'Restart',okay=False)
        model['time']['control']='adaptive_flow'
        rows=run(scheme+'adaptive')
        audit=rows[-1]['terminal_physical_audit']['committed_convective_cfl']
        assert abs(audit['out_max']-.3)<1e-11,audit
        assert abs(audit['directional_max']-.1)<1e-11,audit
        model['time']['control']='fixed'
    print('CFL policy: CN/BE and BE/PISO, raw metrics, admission, adaptive dt, 2-to-4 restart, history drift and witness tamper PASS')
