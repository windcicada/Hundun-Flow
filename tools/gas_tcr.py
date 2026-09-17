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
p.add_argument('--cf-double', type=Path)
p.add_argument('--cf-source', type=Path)
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

def run(path, prefix=(), input_text=data, args=()):
    result = subprocess.run(list(prefix)+[str(path.resolve())]+list(args), input=input_text.encode(),
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
mix_rows = []
for case in range(67):
    density = [1. if case < 3 else rng.randrange(8, 65)/32. for _ in range(64)]
    volume = [1. if case < 3 else rng.randrange(8, 65)/32. for _ in range(64)]
    scalar = []
    for z in range(4):
        for y in range(4):
            for x in range(4):
                scalar.append(.25 if case == 0 else .125*x+.0625*y+.03125*z if case == 1 else
                              (x*x+2*y*y+3*z*z)/64. if case == 2 else rng.randrange(8, 57)/64.)
    mix_rows.append(density+volume+scalar)
mix_data = ''.join(' '.join(format(v, '.17g') for v in row)+'\n' for row in mix_rows)
mix_native = run(a.hundun, shlex.split(a.runner), mix_data, ('mix',))
mix_errors = {}
for precision, path in [('fp64', a.fp64), ('fp32', a.fp32)]:
    reference = run(path, input_text=mix_data, args=('mix',))
    assert len(reference) == len(mix_native) == 8*len(mix_rows)
    maxima = [0., 0.]
    for left, right in zip(mix_native, reference):
        assert len(left) == len(right) == 2
        for col, (x, y) in enumerate(zip(left, right)):
            maxima[col] = max(maxima[col], abs(x-y)/max(1., abs(x), abs(y)))
    mix_errors[precision] = dict(raw_ratio=maxima[0], cphi=maxima[1])
print('mixture_fraction_filter', json.dumps(mix_errors, indent=2))
assert max(mix_errors['fp64'].values()) <= 1e-11, mix_errors
assert max(mix_errors['fp32'].values()) <= 2e-4, mix_errors
cf_filter = None
if a.cf_double:
    assert a.cf_source and a.cf_source.is_file()
    assert sha(a.cf_source) == 'ba15122b54e135ebe1bd15eab7ba351f9c37ff324497bfca64462ee6e1213d9e'
    native_cf = run(a.hundun, shlex.split(a.runner), mix_data, ('cf_mix',))
    reference_cf = run(a.cf_double, input_text=mix_data, args=('cf_mix',))
    assert len(native_cf) == len(reference_cf) == 8*len(mix_rows)
    assert all(len(left) == len(right) == 3 for left,right in zip(native_cf,reference_cf))
    cf_error = max(abs(x-y)/max(1.,abs(x),abs(y)) for left,right in zip(native_cf,reference_cf)
                   for x,y in zip(left,right))
    print('624CF complete Dynamic_Cphi FP64 difference', cf_error)
    assert cf_error <= 1e-11
    cf_filter = dict(cases=len(mix_rows), difference=cf_error, source_sha256=sha(a.cf_source),
        binary_sha256=sha(a.cf_double),
        scope='Complete reference Dynamic_Cphi; equal scalar channels with unit molecular weights isolate moments and bounded products')
# Non-unit molecular weights expose a distinct source-coordinate policy.
# Compare the literal source first, then measure the dimensional correction;
# the latter is a documented model difference rather than an equivalence test.
species_results = {}
for label, path, mode in [('dyn711', a.fp64, 'species'), ('624cf', a.cf_double, 'cf_species')]:
    if not path:
        continue
    literal = run(a.hundun, shlex.split(a.runner), mix_data, (mode,))
    mass = run(a.hundun, shlex.split(a.runner), mix_data, (mode+'_mass',))
    reference = run(path, input_text=mix_data, args=(mode,))
    assert len(reference) == len(literal) == len(mass) == 8*len(mix_rows)
    width = 4 if label == 'dyn711' else 3
    assert all(len(x) == len(y) == len(z) == width for x,y,z in zip(reference,literal,mass))
    error = max(abs(x-y)/max(1.,abs(x),abs(y)) for left,right in zip(literal,reference)
                for x,y in zip(left,right))
    assert error <= 1e-11, (label,error)
    columns = [3] if label == 'dyn711' else [0,1,2]
    differences = [dict(case=i//8, cell=i%8, channel=j, source=left[j], mass_consistent=right[j],
                        absolute_difference=abs(left[j]-right[j]))
                   for i,(left,right) in enumerate(zip(reference,mass)) for j in columns]
    worst = max(differences,key=lambda x:x['absolute_difference'])
    species_results[label] = dict(cases=len(mix_rows), source_equivalent_difference=error,
        changed_coefficients=sum(x['absolute_difference']>1e-11 for x in differences),
        coefficient_count=len(differences), worst=worst,
        molecular_weights=([1.,16.,17.] if label=='dyn711' else [18.,32.,17.]),
        policy='source G2 uses specific-mole gradient; RD2G2 uses mass gradient; mass_consistent uses mass gradients in both')
print('species_coordinate_policy', json.dumps(species_results,indent=2))
evidence = dict(schema='hundun.gas.tcr.v1', calls=len(rows), independent_windows=40*3,
    reference_scope='Complete statistics.F90 with uniform 2x2x2 interior, serial reduction and Cartesian gradient fixture',
    verified=['eight-interval signed accumulation', 'separate ninth evaluation call',
              'kappa lower/upper/weak-rate branches', 'time-scale branches on uniform fields',
              'Cphi update clock and uniform fallback',
              'nonuniform mixture-fraction Cphi filters, priority and smoothing'],
    remaining=['production species-coordinate policy', 'production spatial halos and IBM',
               'production history and scheduling',
               'remixing and gas coupling', 'parallel restart and real-case trajectories'],
    errors=errors, input_sha256=hashlib.sha256(data.encode()).hexdigest(),
    spatial_filter=dict(cases=len(mix_rows), interior_cells=8, errors=mix_errors,
        input_sha256=hashlib.sha256(mix_data.encode()).hexdigest(),
        scope='Complete Dynamic_Cphi mixture-fraction channel; uniform donor units and varying density/volume/scalar; explicit fixture ghost products'),
    cf_filter=cf_filter, species_coordinates=species_results,
    source_sha256=sha(a.source), driver_sha256=sha(Path(__file__).with_suffix('.f90')),
    generator_sha256=sha(Path(__file__)),
    binaries={k:sha(v) for k,v in [('fp32',a.fp32),('fp64',a.fp64),('hundun',a.hundun)]},
    reference_compile=['gfortran -O2 -fcheck=all',
                       'gfortran -O2 -fcheck=all -fdefault-real-8 -fdefault-double-8'],
    identity='dyn711_candidate_kernel; existing cdphyso_dynamic_v1 production semantics preserved')
a.output.write_text(json.dumps(evidence, indent=2)+'\n')
