#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Compare a COAST screen with Hundun monitor records; retain timing scopes."""
import argparse,json,re,hashlib
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--coast',type=Path,required=True)
p.add_argument('--hundun',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
s=a.coast.read_text()
w=[float(x) for x in re.findall(r'wall_step\[s\]\s+([\d.]+)',s)]
starts=[float(x) for x in re.findall(r'Chemistry reactor start time\s*:\s*([\d.]+)',s)]
ends=[float(x) for x in re.findall(r'Chemistry reactor end time\s*:\s*([\d.]+)',s)]
rows=[json.loads(line) for line in a.hundun.read_text().splitlines()]
if not w or not rows or len(starts)!=len(ends):raise ValueError('incomplete timing window')
h=[x['payload']['seconds'] for x in rows]
report=dict(coast=dict(wall_seconds=w,mean_seconds=sum(w)/len(w),
    chemistry_rank0_cpu_seconds=[y-x for x,y in zip(starts,ends)]),
    hundun=dict(steps=[x['step'] for x in rows],dt=[x['payload']['dt'] for x in rows],
    seconds=h,mean_seconds=sum(h)/len(h),records=rows),
    ratio=(sum(h)/len(h))/(sum(w)/len(w)),
    scope='short window; COAST master wall vs Hundun maximum-rank advance; CPU chemistry separate',
    inputs={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [a.coast,a.hundun]})
a.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('hundun','inputs')},indent=2))
