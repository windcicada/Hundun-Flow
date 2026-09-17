#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Compare complete original REAL4 cgsol against the native FP64 IC/PCG probe."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import shlex
import subprocess


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def case(shape, periodic, seed, grounded=False, solid=False):
    rng = random.Random(seed)
    nx, ny, nz = shape
    count = nx * ny * nz
    coords = list(itertools.product(range(nz), range(ny), range(nx)))
    def flat(x, y, z):
        return x + nx * (y + ny * z)
    neighbours = [[None] * 6 for _ in range(count)]
    rows = [[0.] * 6 for _ in range(count)]
    for z, y, x in coords:
        i = flat(x, y, z)
        for axis, length in enumerate(shape):
            at = [x, y, z]
            at[axis] += 1
            if at[axis] == length:
                if not periodic[axis]:
                    continue
                at[axis] = 0
            j = flat(*at)
            neighbours[i][2 * axis + 1] = j
            neighbours[j][2 * axis] = i
            a = (2. ** (2 * axis)) * rng.randint(1, 8) / 16
            if solid and (i == count // 2 or j == count // 2):
                a = 0.
            rows[i][2 * axis + 1] = rows[j][2 * axis] = a
    exact = [(i % 13 - 6) / 8. for i in range(count)]
    data = []
    for i in range(count):
        isolated = solid and i == count // 2
        surplus = (1. if isolated else 0.125 if not grounded or i == 0 else 0.)
        diagonal = sum(rows[i]) + surplus
        rhs = diagonal * exact[i] - sum(a * exact[j] for a, j in zip(rows[i], neighbours[i]) if a)
        volume = 2. ** ((i * 7 + seed) % 9 - 4)
        data.append([volume, diagonal / volume] + [a / volume for a in rows[i]] + [rhs / volume])
    lines = [' '.join(map(str, shape + periodic))]
    lines += [' '.join(format(v, '.17g') for v in row) for row in data]
    return {'shape': shape, 'periodic': periodic, 'seed': seed, 'grounded': grounded,
            'solid': solid, 'rows': data, 'neighbours': neighbours, 'exact': exact,
            'input': '\n'.join(lines) + '\n'}


def invoke(command, cases):
    proc = subprocess.run(command, input=''.join(c['input'] for c in cases),
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                          timeout=180)
    if proc.returncode:
        raise RuntimeError('{} failed: {}\n{}'.format(command, proc.stderr, proc.stdout[-2000:]))
    lines = iter(proc.stdout.splitlines())
    results = []
    for item in cases:
        head = next(lines).split()
        rows = [list(map(float, next(lines).split())) for _ in item['rows']]
        if any(len(row) != 10 or not all(math.isfinite(x) for x in row) for row in rows):
            raise RuntimeError('invalid probe rows')
        results.append({'iterations': int(head[0]), 'reported_residual': float(head[1]), 'rows': rows})
    if next(lines, None) is not None:
        raise RuntimeError('trailing probe output')
    return results


def metrics(item, result, reference=False):
    exact = item['exact']
    xs = [r[0] for r in result['rows']]
    residual = 0.
    matrix = 0.
    rhs_error = 0.
    for i, source in enumerate(item['rows']):
        volume, diagonal = source[:2]
        ax = diagonal * xs[i] - sum(a * xs[j] for a, j in zip(source[2:8], item['neighbours'][i]) if a)
        residual = max(residual, abs(ax - source[8]))
        expected_rhs = source[8] * (1. if reference else volume)
        rhs_error = max(rhs_error, abs(result["rows"][i][9] - expected_rhs) / max(1., abs(expected_rhs)))
        # Original cgsol keeps bpc unscaled and stores its scaled RHS in r.
        expected = [diagonal * volume] + [a * volume for a in source[2:8]]
        for a, b in zip(result['rows'][i][2:9], expected):
            matrix = max(matrix, abs(a-b) / max(1., abs(b)))
    return {'iterations': result['iterations'], 'reported_residual': result['reported_residual'],
            'original_row_residual_max': residual, 'scaled_matrix_difference': matrix, 'rhs_storage_difference': rhs_error,
            'solution_error': max(abs(a-b) for a, b in zip(xs, exact)) / max(1., max(map(abs, exact)))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', required=True)
    parser.add_argument('--hundun', required=True)
    parser.add_argument('--runner', default='')
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    cases = []
    for seed, shape in enumerate([(1, 1, 1), (7, 1, 1), (6, 5, 1), (6, 5, 4), (3, 7, 5)]):
        for periodic in [(0, 0, 0), (1, 0, 0), (1, 1, 1)]:
            cases.append(case(shape, periodic, seed + 11))
    for shape in [(6, 5, 4), (3, 7, 5)]:
        cases.append(case(shape, (1, 0, 1), 41, grounded=True))
        cases.append(case(shape, (1, 0, 1), 43, solid=True))
    reference = invoke([args.reference], cases)
    native = invoke(shlex.split(args.runner) + [args.hundun], cases)
    results = []
    for item, ref, cpp in zip(cases, reference, native):
        one = {key: item[key] for key in ('shape', 'periodic', 'seed', 'grounded', 'solid')}
        one['input_sha256'] = hashlib.sha256(item['input'].encode()).hexdigest()
        one['reference_real4'] = metrics(item, ref, reference=True)
        one['hundun_fp64'] = metrics(item, cpp)
        one['inverse_pivot_relative_difference'] = max(abs(r[1]-c[1])/abs(c[1]) for r, c in zip(ref['rows'], cpp['rows']))
        if one['hundun_fp64']['solution_error'] > 1e-11:
            raise RuntimeError('FP64 solution threshold: {}'.format(one))
        if one['reference_real4']['solution_error'] > 1e-4 or one['inverse_pivot_relative_difference'] > 3e-6:
            raise RuntimeError('REAL4 comparison threshold: {}'.format(one))
        if max(one['reference_real4']['scaled_matrix_difference'], one['hundun_fp64']['scaled_matrix_difference']) != 0:
            raise RuntimeError('dyadic volume scaling differs: {}'.format(one))
        if max(one['reference_real4']['rhs_storage_difference'], one['hundun_fp64']['rhs_storage_difference']) != 0:
            raise RuntimeError('RHS storage differs: {}'.format(one))
        results.append(one)
    output = {'schema': 'hundun-gas-iccg-reference-v1', 'scope': 'complete unchanged serial cgsol, shared dyadic SPD matrices; native PCG true residual and IC(0)',
              'reference_precision': 'REAL(kind=4), including explicit work pointers',
              'native_precision': 'FP64', 'case_count': len(cases),
              'controls': {'reference_original_row_max_tolerance': 2e-6, 'reference_rnorm': 1.,
                           'native_scaled_l2_atol': 1e-12, 'native_scaled_l2_rtol': 1e-13,
                           'maximum_iterations': 2000, 'purpose': 'correctness; iteration counts are not performance-equivalent stopping criteria'},
              'source_sha256': digest(args.source),
              'driver_sha256': digest(Path(__file__).with_suffix('.f90')),
              'generator_sha256': digest(__file__),
              'reference_build': 'mpifort -O2 -cpp -Jcheck/gi tools/gas_iccg.f90 check/gi/cgsol.F90 -o check/gi/ref',
              'reference_rhs_storage': 'coef(bpc) retains original RHS; initial volume-scaled r is subsequently updated by CG',
              'native_rhs_storage': 'scaled row RHS', 'reference_binary_sha256': digest(args.reference),
              'native_binary_sha256': digest(args.hundun), 'results': results}
    Path(args.output).write_text(json.dumps(output, indent=2) + '\n')
    print('ICCG {} cases; FP64 error {:.3e}; REAL4 error {:.3e}; pivot difference {:.3e}'.format(
        len(cases), max(r['hundun_fp64']['solution_error'] for r in results),
        max(r['reference_real4']['solution_error'] for r in results),
        max(r['inverse_pivot_relative_difference'] for r in results)))


if __name__ == '__main__':
    main()
