#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Compare complete cmod/step to the native CN row and independent BE/pressure rows."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import shlex
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--float', dest='fp32', type=Path, required=True)
p.add_argument('--double', dest='fp64', type=Path, required=True)
p.add_argument('--hundun', type=Path, required=True)
p.add_argument('--runner', default='')
p.add_argument('--cmod', type=Path, required=True)
p.add_argument('--step', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
rng = random.Random(711)
# storage, dt, old centre/neighbours, spatial diagonal/neighbours/RHS, drho/dt.
rows = [[2., .1] + [3.]*7 + [6.] + [1.]*6 + [4., .2]]
for i in range(256):
    neighbours = [rng.uniform(0, 100) for _ in range(6)]
    rows.append([10**rng.uniform(-2, 2), 10**rng.uniform(-5, 0)] +
                [rng.uniform(-100, 100) for _ in range(7)] +
                [sum(neighbours)+rng.uniform(0, 5)] + neighbours +
                [rng.uniform(-100, 100), rng.uniform(-2, 2)])
def text(inputs):
    return ''.join(' '.join(format(v, '.17g') for v in r)+'\n' for r in inputs)
def run(path, data, prefix=()):
    result = subprocess.run(list(prefix)+[str(path.resolve())], input=data.encode(),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return [list(map(float, line.split())) for line in result.stdout.decode().splitlines()]
reference_input = text([[mode]+row for row in rows for mode in (1, 2, 3)])
f4 = run(a.fp32, reference_input)
f8 = run(a.fp64, reference_input)
hf = run(a.hundun, text([row[:-1] for row in rows]), shlex.split(a.runner))
assert len(f4) == len(f8) == 3*len(hf) == 3*len(rows)
errors = {k: 0. for k in ('cn_fp64', 'cn_fp32', 'be_fp64', 'be_fp32', 'pressure_fp64', 'pressure_fp32')}
for i, row in enumerate(rows):
    rho, dt = row[:2]
    old = row[2:9]
    spatial = row[9:17]
    be = list(spatial)
    be[0] += rho/dt
    be[-1] += rho/dt*old[0]
    pressure = list(spatial)
    pressure[-1] -= row[-1]
    for mode, name, expected in ((0, 'cn', hf[i]), (1, 'be', be), (2, 'pressure', pressure)):
        # RHS cancellation is normalized by the magnitude of its input terms.
        # Matrix coefficients retain their individual relative scale.
        scale = [max(1e-30, abs(v)) for v in expected]
        scale[-1] = max(1e-30, abs(spatial[-1])+abs(rho/dt*old[0])+
                        sum(abs(x*y) for x, y in zip(spatial[:7], old))+abs(row[-1]))
        for precision, result in (('fp64', f8[3*i+mode]), ('fp32', f4[3*i+mode])):
            error = max(abs(x-y)/s for x, y, s in zip(expected, result, scale))
            key = name+'_'+precision
            errors[key] = max(errors[key], error)
assert all(v < (1e-11 if k.endswith('fp64') else 8e-7) for k, v in errors.items()), errors
# Constant old fields remove a conservative spatial row; the source remains full.
assert hf[0] == [23.] + [.5]*6 + [64.], hf[0]
# Same accepted history and identical frozen row produce the same repeated solve.
repeat = run(a.hundun, text([rows[1][:-1], rows[1][:-1]]), shlex.split(a.runner))
assert repeat[0] == repeat[1]
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
evidence = dict(schema='hundun.gas.time.v1', samples=len(rows), modes=['CN', 'BE', 'pressure'],
                errors=errors, input_sha256=hashlib.sha256(reference_input.encode()).hexdigest(),
                source_sha256={'cmod': sha(a.cmod), 'step': sha(a.step)},
                binaries={k: sha(v) for k, v in [('fp32', a.fp32), ('fp64', a.fp64), ('hundun', a.hundun)]},
                coverage='complete cmod/step; native CN row; analytic BE/pressure; isolated frozen row only',
                result='pass')
a.output.write_text(json.dumps(evidence, indent=2)+'\n')
print(json.dumps(evidence))
