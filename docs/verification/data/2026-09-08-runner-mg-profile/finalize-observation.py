#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Read-only cross-run checks and a receipt for the already completed MG pilot."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

audit = Path(__file__).resolve().parent
base = audit.parent
new = base / 'pilot-mg-9500-9510-20260908'
old = base / 'pilot-criterion-9500-9510-20260908'
checkout = Path('/home/wyf/code_dev/.worktrees/hundun-flow-mg-runner-accept-20260908')
output = audit / 'PILOT_ACCEPTED.json'
assert not output.exists(), 'refusing to replace an existing receipt'

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()

def generation(root):
    name = (root / 'Restart/current').read_text().strip()
    assert name.startswith('generation-9510-') and '/' not in name
    return root / 'Restart' / name

def equal_group(left, right, pattern, expected):
    names = sorted(p.name for p in left.glob(pattern))
    assert names == sorted(p.name for p in right.glob(pattern))
    assert len(names) == expected
    manifest = hashlib.sha256()
    for name in names:
        a, b = left / name, right / name
        assert a.stat().st_size == b.stat().st_size, name
        with a.open('rb') as x, b.open('rb') as y:
            while True:
                chunk = x.read(1048576)
                assert chunk == y.read(1048576), name
                if not chunk:
                    break
        manifest.update((sha(a) + '  ' + name + '\n').encode('ascii'))
    return {'files': len(names), 'byte_equal': True,
            'sorted_sha256_manifest_digest': manifest.hexdigest()}

for checks in ('SOURCE.sha256', 'FROZEN.sha256'):
    subprocess.run(['sha256sum', '--check', '--status', checks], cwd=str(audit), check=True)
local = json.loads((audit / 'LOCAL_ACCEPTED.json').read_text())
observation = json.loads((audit / 'observation.json').read_text())
physical = json.loads((audit / 'pilot-output-verified.json').read_text())
assert observation['complete'] and physical['complete']
assert physical['steps'] == 10 and physical['ranks'] == 128
ng, og = generation(new), generation(old)
equivalence = {
    'comparison_root': str(old),
    'rank_checkpoint': equal_group(ng, og, 'rank-*.bin', 128),
    'manifest': equal_group(ng, og, 'manifest.bin', 1),
    'visit': equal_group(new / 'Visit', old / 'Visit', '*.vtr', 128),
    'accumulator': equal_group(new, old, 'step-*.accumulator', 1),
    'statistics': equal_group(new, old, 'step-*.statistics.json', 1),
}
def completion(root):
    path = root / 'step-00000000000000009510.complete'
    lines = path.read_text().splitlines()
    record = [line for line in lines if line.startswith('restart_generation ')]
    assert record == ['restart_generation ' + generation(root).name]
    return [line for line in lines if not line.startswith('restart_generation ')]
assert completion(new) == completion(old)
equivalence['completion_difference'] = 'Only restart_generation differs; each names its own verified generation.'

manifest = base / 'long-thermodynamic-35000-20260907/Restart/generation-9500-139424060480232/manifest.bin'
command = [sys.executable, str(checkout / 'tools/v04_evidence_validate.py'),
           'runtime', str(new / 'evidence.jsonl'), '--run-start-manifest', str(manifest)]
runtime = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True, timeout=60)
assert runtime.returncode == 0, runtime.stderr
test_log = checkout / 'build-accept/Testing/Temporary/LastTest.log'
assert test_log.read_text().count('\nTest Passed.\n') == 13
assert '\nTest Failed.' not in test_log.read_text()
launch_log = new.with_suffix('.log')
assert 'COMPLETED steps=10 final_step=9510 samples=0' in launch_log.read_text()
receipt = {
    'schema': 'HUNDUN_MG_OBSERVATION_PILOT_ACCEPTANCE_V1',
    'scope': 'Observation-only acceptance, not an algorithm speedup or COAST replacement acceptance.',
    'source_commit': local['source_commit'], 'source_tree': local['source_tree'],
    'runner_sha256': local['runner_sha256'],
    'runner_build_manifest_sha256': local['runner_build_manifest_sha256'],
    'retained_method_signature': local['retained_method_signature'],
    'source_generation': local['source_generation'], 'output_generation': ng.name,
    'window': local['pilot'], 'process_wall_seconds': 110.86,
    'clean_release_tests': local['clean_release_tests'],
    'physical_audit': {k: physical[k] for k in (
        'steps', 'ranks', 'bdf1_steps', 'bdf2_steps', 'terminal_maxima', 'fluid_ranges')},
    'payload_equivalence': equivalence,
    'runtime_validation': {'command': command, 'exit_code': runtime.returncode,
                           'stdout': runtime.stdout, 'stderr': runtime.stderr},
    'source_files_unchanged': local['source_integrity_files'],
    'artifacts_sha256': {str(p.relative_to(base)) if base in p.parents else str(p): sha(p)
                       for p in [Path(__file__), audit / 'summarize-cost.py',
                                 audit / 'cost-summary.json', audit / 'observation.json',
                                 audit / 'loops.jsonl', audit / 'pilot-output-verified.json',
                                 audit / 'SOURCE.sha256', audit / 'FROZEN.sha256',
                                 test_log, launch_log]},
    'limitations': ['Single ten-step observation window; no timing A/B or speedup claim.',
                    'Solid regional extrema checked; no new source-to-final solid-only payload comparison.',
                    'Original long-run processes remain paused beyond this pilot endpoint; do not replace them.',
                    'Stage 6/two-phase integration remains deferred.'],
}
with output.open('x') as stream:
    json.dump(receipt, stream, indent=2, ensure_ascii=True)
    stream.write('\n')
print('PASS: 128 checkpoint ranks, 128 Visit files, manifest, accumulator and statistics byte equal; runtime exit 0.')
print(str(output))
