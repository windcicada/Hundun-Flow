#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary, mpi=map(lambda p:str(Path(p).resolve()),sys.argv[1:3])
def call(n,args,cwd,ok=True):
    r=subprocess.run([mpi,'--oversubscribe','--bind-to','none','-n',str(n),binary,*map(str,args)],
                     cwd=str(cwd),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                     universal_newlines=True,timeout=45)
    assert (r.returncode==0)==ok,r.stdout[-6000:]
    return r.stdout
with tempfile.TemporaryDirectory(prefix='hf-run-') as tmp:
    root=Path(tmp);case=root/'c'
    call(1,['init-case','--output',case],root)
    config=dict(mode='new',steps=2,output='../r',monitor_interval=1,
                output_interval=0,restart_interval=1,max_dt=1e-5)
    (case/'run.json').write_text(json.dumps(config))
    text=call(2,[],case)
    assert 'COMPLETED steps=2 ' in text
    rows=[json.loads(s) for s in (root/'r/monitor.jsonl').read_text().splitlines()]
    assert [r['step'] for r in rows]==[1,2]
    assert not (root/'r/Visit').exists()
    text=call(2,['run',case,'--config',case/'run.json','--steps',1,'--output',root/'s'],root)
    assert 'COMPLETED steps=1 ' in text
    config.update(mode='restart',restart='../r/Restart',output='../rr')
    del config['steps'];config['end_time']=2.5e-5
    (case/'run.json').write_text(json.dumps(config))
    assert 'COMPLETED steps=3 ' in call(4,[],case)
    bad=['{"mode":"new","mode":"restart","steps":1,"output":"../bad"}',
         '{"mode":"new","steps":1,"output":"../bad","unknown":1}',
         '{"mode":"restart","steps":1,"output":"../bad"}',
         '{"mode":"new","steps":1,"output":"../bad","restart":"../r/Restart"}']
    for s in bad:
        (case/'run.json').write_text(s);call(2,[],case,False)
    assert not (root/'bad').exists()
print('run.json: default/explicit, MPI, CLI overrides, end time, restart and strict schema PASS')
