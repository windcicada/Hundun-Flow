#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Native and statistics entry points share CFL stepping and exact continuation."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--native', type=Path, required=True)
p.add_argument('--runner', type=Path, required=True)
p.add_argument('--mpi', required=True)
a = p.parse_args()
repo = Path(__file__).resolve().parents[4]
with tempfile.TemporaryDirectory(prefix='hf-cfl-') as tmp:
    root = Path(tmp)
    case = root / 'case'
    case.mkdir()
    m = json.loads((repo / 'examples/minimal/case.json').read_text())
    m['mesh'].update(domain={'lower': [-2, -2, -2], 'upper': [2, 2, 2]},
        exact_cells=[16, 16, 16], minimum_spacing=[.25] * 3,
        immersed_boundary={'stl_file': 'cyl.stl', 'fluid_side': 'outside',
                           'reconstruction_policy': 'adaptive_order'})
    m['mesh']['limits'].update(max_global_cells=4096, max_memory_bytes_per_rank=1073741824)
    m['time'].update(initial_dt=1., maximum_dt=1., convective_cfl=.002)
    m['solver']['cold_stopping'] = dict(reference_time=.002, momentum=1e-4, enthalpy=1e-4, species=1e-4)
    (case / 'case.json').write_text(json.dumps(m))
    shutil.copyfile(repo / 'examples/minimal/thermophysics.d', case / 'thermophysics.d')
    shutil.copyfile(repo / 'versions/v0.4/tests/data/cylinder_ascii.stl', case / 'cyl.stl')
    spec = root / 'stat.d'
    spec.write_text('HUNDUN_V04_LITERATURE_STATISTICS_V1\ndevelopment_steps 0\n'
        'collection_end_step 10\ncheckpoint_interval 2\nrho_ref 1\nu_ref 1\n'
        'diameter 1\nspan 4\ncylinder_center_x 0\nstation_x_over_d 1.5\nend\n')

    def run(binary, args, name, restart=None):
        result = subprocess.run([a.mpi, '-n', '2', str(binary.resolve()), *map(str, args)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=90)
        assert result.returncode == 0, result.stdout[-12000:]
        path = root / name / 'evidence.jsonl'
        command = [sys.executable, str(repo / 'tools/v04_evidence_validate.py'), 'runtime', str(path)]
        if restart:
            generation = (restart / 'current').read_text().strip()
            command += ['--run-start-manifest', str(restart / generation / 'manifest.bin')]
        checked = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        assert checked.returncode == 0, checked.stdout
        return [json.loads(line) for line in path.read_text().splitlines()]

    native = run(a.native, ['run', case, '--output', root / 'native', '--steps', 4,
        '--output-interval', 0, '--restart-interval', 0, '--diagnostics-interval', 0], 'native')
    common = ['--case-root', case, '--spec', spec, '--steps', 2, '--visit-interval', 0]
    stats = run(a.runner, [*common, '--run-root', root / 'stats'], 'stats')
    resumed = run(a.runner, [*common, '--run-root', root / 'resume',
        '--restart-root', root / 'stats/Restart'], 'resume', root / 'stats/Restart')
    for expected, actual in zip(native, stats + resumed):
        assert expected['step'] == actual['step']
        assert expected['time'] == actual['time'], (expected['time'], actual['time'])
        assert actual['coupling'] == 'CN_BE'
        assert expected['retry'] == actual['retry']
        assert actual['cold']['solid_velocity_max'] == 0
    assert len(native) == 4 and len(stats + resumed) == 4
    assert native[0]['time'] < .001  # Initial dt=1 is actually reduced by the field CFL.
    assert native[1]['time'] - native[0]['time'] != native[0]['time']
    print('PASS shared CFL, CN/BE, IBM, variable dt and exact statistics continuation')
