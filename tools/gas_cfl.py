#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Compare the original complete courant routine with Hundun's cell operator."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import shlex
import subprocess

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--float',dest='fp32',required=True,type=Path)
p.add_argument('--double',dest='fp64',required=True,type=Path)
p.add_argument('--hundun',required=True,type=Path)
p.add_argument('--runner',default='')
p.add_argument('--source',required=True,type=Path)
p.add_argument('--output',required=True,type=Path)
a=p.parse_args()
rng=random.Random(711)
rows=[[1,1,.3,1,1,1,1,1,1],[2,3,.25,-1,2,-3,4,-5,6],[1,1,.1,0,0,0,0,0,0]]
for i in range(256):
    rows.append([10**rng.uniform(-2,2),10**rng.uniform(-6,0),10**rng.uniform(-7,-3)]+
                [rng.uniform(-10,10) for j in range(6)])
text=''.join(' '.join(format(v,'.17g') for v in row)+'\n' for row in rows)
def run(path,prefix=()):
    r=subprocess.run(list(prefix)+[str(path.resolve())],input=text.encode(),stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,check=True)
    return [list(map(float,line.split())) for line in r.stdout.decode().splitlines()]
f4=run(a.fp32);f8=run(a.fp64);hf=run(a.hundun,shlex.split(a.runner))
assert len(f4)==len(f8)==len(hf)==len(rows)
errors={'fp32':0.,'fp64':0.,'outgoing':0.,'absolute':0.}
for row,low,high,got in zip(rows,f4,f8,hf):
    for name,ref in [('fp32',low[0]),('fp64',high[0])]:
        errors[name]=max(errors[name],abs(got[0]-ref)/max(1e-30,abs(got[0]),abs(ref)))
    rho,vol,dt,*f=row
    out=sum(max(v*(-1 if i%2==0 else 1),0) for i,v in enumerate(f))*dt/(rho*vol)
    absolute=sum(map(abs,f))*.5*dt/(rho*vol)
    for name,x,y in [('outgoing',got[1],out),('absolute',got[2],absolute)]:
        errors[name]=max(errors[name],abs(x-y)/max(1e-30,abs(x),abs(y)))
assert errors['fp64']<1e-11 and errors['fp32']<5e-7,errors
assert errors['outgoing']<1e-11 and errors['absolute']<1e-11
assert abs(hf[0][0]-.3)<1e-15 and abs(hf[0][1]-.9)<1e-15
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
evidence=dict(schema='hundun.gas.cfl.v1',samples=len(rows),errors=errors,
              input_sha256=hashlib.sha256(text.encode()).hexdigest(),
              source_sha256=sha(a.source),binaries={k:sha(v) for k,v in
                  [('fp32',a.fp32),('fp64',a.fp64),('hundun',a.hundun)]},
              coverage='complete courant.F90; serial identity collective; local operator only',
              result='pass')
a.output.write_text(json.dumps(evidence,indent=2)+'\n')
print(json.dumps(evidence))
