#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Real CLI: detail numbers alone cannot identify the failure category."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--archive-binary', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = args.output or Path(tempfile.mkdtemp(prefix='hundun-failure-hint-')) / 'audit'
    root.mkdir(parents=True, exist_ok=False)
    case = root / 'case'
    case.mkdir()
    data = Path(__file__).resolve().parents[1] / 'data'
    model = json.loads((data / 'case_minimal_valid.json').read_text())
    model['boundaries']['y_min'].update(flow_kind='no_slip_wall',
        thermal_kind='isothermal_wall', temperature=400.0, velocity=[0.0, 0.0, 0.0])
    (case / 'case.json').write_text(json.dumps(model))
    shutil.copyfile(str(data / 'thermophysics.d'), str(case / 'thermophysics.d'))

    def run(binary, name, extra):
        command = [str(binary.resolve()), 'run', str(case.resolve()),
            '--output', str((root / name).resolve()), '--steps', '1',
            '--output-interval', '0', '--restart-interval', '0'] + extra
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                universal_newlines=True, timeout=30)
        (root / (name + '.log')).write_text(result.stdout)
        assert result.returncode != 0, result.stdout
        assert 'COMPLETED steps=' not in result.stdout, result.stdout
        assert 'committed_step=0' in result.stdout and 'attempts=0' in result.stdout, result.stdout
        return result.stdout

    output = run(args.binary, 'ambiguous', [])
    assert 'termination_phase=initialize' in output and 'detail=10505' in output, output
    assert 'conflicting boundary-derived initial state' in output, output
    print('PASS ambiguous initial state keeps the explicit-state hint', flush=True)
    if args.archive_binary:
        # The second executable must actually be built without Git metadata.
        # No product function or identity check is mocked/bypassed here.
        output = run(args.archive_binary, 'archive-identity',
                     ['--initial-state', '101325,300,0.1,0,0'])
        assert 'termination_phase=runtime_identity' in output and 'detail=10505' in output, output
        assert 'conflicting boundary-derived initial state' not in output, output
        print('PASS unavailable archive identity is not mislabeled as initial-state conflict', flush=True)


if __name__ == '__main__':
    main()
