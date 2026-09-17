#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Compare dyn711 window/branches with the complete frozen statistics routine."""
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
p.add_argument('--source', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
assert sha(a.source) == '31208e34befcabdf9fdfcfb744cb8f53b25b4d5db94b329082b6f49acf9d24a3'
rng = random.Random(711)
rows = []
for case in range(40):
    eta = (.25, .5, .75, .03125, .96875)[case % 5]
    for call in range(27):
        # Variable dyadic intervals isolate algebra from parser precision.
        dt = 2.**-rng.randrange(8, 16)
        psr = [2.**rng.randrange(-5, 5) for _ in range(6)]
        pdf = [v*r for v, r in zip(psr, (.5, 1., 2., -2., 0., 1.))]
        psr[3] *= -1  # signed same-direction destruction vs opposite signs
        psr[5] = 0   # weak reference rate: defined unit-ratio branch
        if call % 9 == 8:
            pdf = [1e12]*6  # evaluation-call rates lie outside the window
        epsilon = 2.**(32 if case % 2 else -8)
        velocity = 2.**(-case % 3)
        composition = [.125, .25, .25, .0625, .125, 0.]
        rows.append([int(call == 0), dt, eta, epsilon, velocity]+pdf+psr+composition)
data = ''.join(' '.join(format(v, '.17g') for v in row)+'\n' for row in rows)

def run(path, prefix=()):
    result = subprocess.run(list(prefix)+[str(path.resolve())], input=data.encode(),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return [list(map(float, line.split())) for line in result.stdout.decode().splitlines()]

native = run(a.hundun, shlex.split(a.runner))
results = {'fp64': run(a.fp64), 'fp32': run(a.fp32)}
errors = {}
for precision, reference in results.items():
    assert len(native) == len(reference) == len(rows)
    worst = (0., 0, 0)
    for row, (left, right) in enumerate(zip(native, reference)):
        assert len(left) == len(right) == 30
        assert left[:4] == right[:4], (precision, row, left[:4], right[:4])
        for col, (x, y) in enumerate(zip(left[4:], right[4:]), 4):
            # Relative for nonzero values; unit floor for near-zero roots.
            error = abs(x-y)/max(1., abs(x), abs(y))
            if error > worst[0]:
                worst = (error, row, col)
    errors[precision] = dict(max_normalized_difference=worst[0], row=worst[1], column=worst[2])
print(json.dumps(errors, indent=2))
assert errors['fp64']['max_normalized_difference'] <= 1e-11, errors
assert errors['fp32']['max_normalized_difference'] <= 5e-7, errors
evidence = dict(schema='hundun.gas.tcr.v1', calls=len(rows), independent_windows=40*3,
    reference_scope='Complete statistics.F90 with uniform 2x2x2 interior, serial reduction and Cartesian gradient fixture',
    verified=['eight-interval signed accumulation', 'separate ninth evaluation call',
              'kappa lower/upper/weak-rate branches', 'time-scale branches on uniform fields',
              'Cphi update clock and uniform fallback'],
    remaining=['nonuniform dynamic Cphi spatial filters', 'production history and scheduling',
               'remixing and gas coupling', 'parallel restart and real-case trajectories'],
    errors=errors, input_sha256=hashlib.sha256(data.encode()).hexdigest(),
    source_sha256=sha(a.source), driver_sha256=sha(Path(__file__).with_suffix('.f90')),
    generator_sha256=sha(Path(__file__)),
    binaries={k:sha(v) for k,v in [('fp32',a.fp32),('fp64',a.fp64),('hundun',a.hundun)]},
    reference_compile=['gfortran -O2 -fcheck=all',
                       'gfortran -O2 -fcheck=all -fdefault-real-8 -fdefault-double-8'],
    identity='dyn711_candidate_kernel; existing cdphyso_dynamic_v1 production semantics preserved')
a.output.write_text(json.dumps(evidence, indent=2)+'\n')
