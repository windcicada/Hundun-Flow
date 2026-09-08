#!/usr/bin/env python3
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Resume only the known frozen job; always pause it again at the handoff."""
import datetime
import fcntl
import pathlib
import subprocess
import time

BASE = pathlib.Path('/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52')
SOURCE = BASE / 'long-thermodynamic-35000-20260907'
SERVICE = 'hundun-re3900-thermodynamic-long-20260907.service'
TARGET = SOURCE / 'step-00000000000000010000.complete'


def command(*args):
    return subprocess.check_output(args, universal_newlines=True).strip()


def signal_job(signal):
    subprocess.check_call(['systemctl', '--user', 'kill', '--kill-who=all', '--signal=' + signal, SERVICE])


def event(message):
    print(datetime.datetime.now().isoformat() + ' ' + message, flush=True)


with (BASE / '.hundun-mpi-maintenance.lock').open('a') as lock:
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not TARGET.exists(), 'target already exists; refuse ambiguous handoff'
    assert command('systemctl', '--user', 'show', SERVICE, '--property=MainPID', '--value') == '183061'
    rows = command('ps', '-eo', 'pid=,stat=,comm=').splitlines()
    ranks = []
    for row in rows:
        pid, state, name = row.split()
        assert name not in ('ctest', 'ninja', 'make', 'clang', 'clang-15', 'cc1plus'), 'compiler/test active'
        if name in ('mpirun', 'mpiexec', 'v04_thin_domain', 'hundun', 'coast', 'COAST'):
            assert state.startswith('T'), 'competing runnable process ' + pid
            assert SERVICE in pathlib.Path('/proc/' + pid + '/cgroup').read_text(), 'unowned process ' + pid
            if name == 'v04_thin_domain':
                ranks.append(pid)
    assert len(ranks) == 128
    assert (SOURCE / 'Restart/current').read_text().strip() == 'generation-9500-139424060480232'
    event('RESUME original 128 ranks from in-memory step 9951; target durable 10000; timeout 1800s')
    try:
        signal_job('SIGCONT')
        deadline = time.monotonic() + 1800
        while not TARGET.exists():
            assert time.monotonic() < deadline, 'bridge timeout'
            assert command('systemctl', '--user', 'is-active', SERVICE) == 'active', 'old job stopped before handoff'
            time.sleep(2)
        lines = TARGET.read_text().splitlines()
        assert lines[0] == 'HUNDUN_V04_THIN_DOMAIN_CHECKPOINT_V1' and lines[-1] == 'end'
        values = dict(line.split(' ', 1) for line in lines[1:-1])
        assert values['step'] == '10000'
        generation = (SOURCE / 'Restart/current').read_text().strip()
        assert values['restart_generation'] == generation and generation.startswith('generation-10000-')
        assert len(list((SOURCE / 'Restart' / generation).glob('rank-*.bin'))) == 128
        for key in ('statistics', 'accumulator'):
            assert (SOURCE / values[key]).stat().st_size > 0
        event('DURABLE_HANDOFF ' + generation)
    finally:
        signal_job('SIGSTOP')
        event('PAUSE original service; physical state is not rolled back; old output is a stopped prefix, not COMPLETED')
    time.sleep(1)
    assert all(command('ps', '-p', pid, '-o', 'stat=').startswith('T') for pid in ranks)
    event('VERIFIED all 128 original ranks stopped')
