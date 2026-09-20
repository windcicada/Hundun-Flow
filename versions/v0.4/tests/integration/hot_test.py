#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Bounded IEM/JL4 execution and native Restart; physical balances are reported."""
import json,math,shutil,subprocess,sys
from pathlib import Path
binary,fixture,mpi,work,validator,compare=map(Path,sys.argv[1:])
if work.exists():shutil.rmtree(work)
shutil.copytree(fixture,work)
model=json.loads((work/'case.json').read_text())
model['reaction']['ensemble']['tcr']={'mode':'off'}
(work/'case.json').write_text(json.dumps(model,indent=2)+'\n')
def call(args,log):
    with log.open('w') as output:
        result=subprocess.run(list(map(str,args)),stdout=output,stderr=output)
    assert result.returncode==0,log.read_text()[-8000:]
def run(name,ranks,restart=None):
    out=work.with_name(work.name+name)
    if out.exists():shutil.rmtree(out)
    args=[mpi,'--oversubscribe','--bind-to','none','-n',ranks,binary,'run',work,
          '--output',out,'--steps',2,'--output-interval',0,'--restart-interval',2,
          '--diagnostics-interval',1]
    args+=['--restart',restart/'Restart'] if restart else ['--initial-state','100000,1100,10,0,0,0.05,0.2,0.01,0.005,0.05,0.001']
    call(args,work/(name+'.log'))
    audit=[sys.executable,validator,'runtime',out/'evidence.jsonl']
    if restart:audit+=['--run-start-manifest',restart/'Restart'/(restart/'Restart/current').read_text().strip()/'manifest.bin']
    call(audit,work/(name+'-audit.log'))
    records=[json.loads(x)['payload'] for x in (out/'diagnostics.jsonl').read_text().splitlines()]
    for row in records:
        assert all(math.isfinite(row[k]) for k in ['mass_kg','internal_energy_J','mass_balance_defect_kg_s','total_energy_balance_defect_W'])
        assert row['mass_kg']>0
    return out
start=run('1',1)
one=run('r',1,start)
two=run('2',2,start)
four=run('4',4,start)
for name,other in [('2',two),('4',four)]:
    call([mpi,'--oversubscribe','--bind-to','none','-n',2,compare,work,one/'Restart',other/'Restart','--reference-linear'],work/(name+'-compare.log'))
print('IEM/JL4 batch run and 1/2/4-rank Restart passed; balance observations retained')
