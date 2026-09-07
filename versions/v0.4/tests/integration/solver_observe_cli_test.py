#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
"""Exercise the public observer CLI, including old/new loop CSV contracts."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

SOURCE = 'source_meta_sha256'


def write(path, records):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def check_linear_criterion(script, checks, root, base, total):
    root.mkdir()
    (root / "conservation.csv").write_text("step,rank\n7,0\n")
    criterion = dict(linear_criterion_valid=1, linear_rhs_norm=100.0,
                     linear_atol=0.0001, linear_rtol=0.001,
                     linear_residual_limit=0.1)
    observed = dict(base, linear_initial=0.5, **criterion)

    def run(rows, schema=4, accepted=True, reason=None):
        metadata = ("HUNDUN_V04_THIN_DOMAIN_RUN_V1\nstarting_step 6\n"
                    "requested_steps 1\nexpected_ranks {}\nobservation_schema {}\n"
                    "end\n").format(len(rows), schema)
        (root / "RUN.meta").write_text(metadata)
        identity = hashlib.sha256(metadata.encode()).hexdigest()
        for path in root.glob("solver-rank-*.csv"):
            path.unlink()
        for rank, row in enumerate(rows):
            write(root / "solver-rank-{}.csv".format(rank),
                  [dict(row, rank=rank, source_meta_sha256=identity)])
        write(root / "performance.csv", [dict(total, rank=rank,
              pressure_calls=row["invoked"], A_apply_ns=row["A_ns"],
              M_apply_ns=row["M_ns"], source_meta_sha256=identity)
              for rank, row in enumerate(rows)])
        result = subprocess.run([sys.executable, script, str(root)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=10)
        assert (result.returncode == 0) == accepted, result.stderr or result.stdout
        if reason:
            assert reason in result.stderr, result.stderr
        checks.append("criterion accepted" if accepted else "criterion rejected")
        return json.loads(result.stdout) if accepted else None

    # The public threshold uses the actual RHS, not the warm-start residual.
    result = run([observed])
    loop = result["loops"][0]
    assert result["schema"] == "HUNDUN_LOOP_OBSERVATION_V4" and result["complete"]
    assert loop["linear_criterion_available"]
    assert loop["linear_criterion"] == dict(rhs_norm=100.0, atol=0.0001,
                                            rtol=0.001, residual_limit=0.1)
    assert loop["linear_contraction"] == [0.02]
    sidecar = root / "criterion-details.jsonl"
    streamed = subprocess.run([sys.executable, script, str(root),
        "--details-output", str(sidecar)], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, universal_newlines=True, timeout=10)
    assert streamed.returncode == 0, streamed.stderr
    streamed = json.loads(streamed.stdout)
    assert streamed["complete"] and streamed["loops"] == []
    assert [json.loads(line) for line in sidecar.read_text().splitlines()] == result["loops"]
    assert streamed["details"]["sha256"] == hashlib.sha256(sidecar.read_bytes()).hexdigest()
    checks.append("criterion streamed details match")

    # Old frozen observations stay complete without inventing a criterion.
    old = run([base], schema=3)
    assert old["schema"] == "HUNDUN_LOOP_OBSERVATION_V3" and old["complete"]
    assert not old["loops"][0]["linear_criterion_available"]
    assert old["loops"][0]["linear_criterion"] is None
    extended = run([observed], schema=3)
    assert extended["loops"][0]["linear_criterion"] == loop["linear_criterion"]
    run([base], accepted=False, reason="criterion")
    for key in criterion:
        incomplete = {k: v for k, v in observed.items() if k != key}
        run([incomplete], accepted=False, reason="criterion")
        run([incomplete], schema=3, accepted=False, reason="criterion")

    unavailable = dict(base, **dict.fromkeys(criterion, 0))
    for invoked in (0, 1):
        missing = run([dict(unavailable, invoked=invoked)])
        assert missing["complete"]
        assert not missing["loops"][0]["linear_criterion_available"]
        assert missing["loops"][0]["linear_criterion"] is None
    for key in criterion:
        if key != "linear_criterion_valid":
            run([dict(unavailable, **{key: criterion[key]})], accepted=False, reason="criterion")
    for key in criterion:
        for value in ("nan", "inf", -1):
            run([dict(observed, **{key: value})], accepted=False)
    for key, value in (("linear_criterion_valid", 2),
                       ("linear_residual_limit", 0), ("linear_residual_limit", 0.10001),
                       ("linear_residual_limit", 0.0005)):
        run([dict(observed, **{key: value})], accepted=False, reason="criterion")
    # Match the public solver controls, which allow either tolerance to be
    # zero and do not cap rtol. A zero RHS can then have a zero limit.
    for edge, expected in (
            (dict(observed, linear_atol=1, linear_rtol=0, linear_residual_limit=1),
             dict(rhs_norm=100.0, atol=1.0, rtol=0.0, residual_limit=1.0)),
            (dict(observed, linear_atol=0),
             dict(rhs_norm=100.0, atol=0.0, rtol=0.001, residual_limit=0.1)),
            (dict(observed, linear_atol=0, linear_rhs_norm=0, linear_residual_limit=0),
             dict(rhs_norm=0.0, atol=0.0, rtol=0.001, residual_limit=0.0)),
            (dict(observed, linear_rtol=1, linear_residual_limit=100),
             dict(rhs_norm=100.0, atol=0.0001, rtol=1.0, residual_limit=100.0)),
            (dict(observed, linear_rtol=2, linear_residual_limit=200),
             dict(rhs_norm=100.0, atol=0.0001, rtol=2.0, residual_limit=200.0))):
        valid = run([edge])
        assert valid["loops"][0]["linear_criterion_available"]
        assert valid["loops"][0]["linear_criterion"] == expected
    run([dict(unavailable, linear_criterion_valid=1)], accepted=False, reason="criterion")
    run([dict(observed, linear_rhs_norm=0, linear_atol=0,
              linear_residual_limit=5e-324)], accepted=False, reason="criterion")
    run([dict(observed, linear_rhs_norm=0, linear_atol=5e-324,
              linear_residual_limit=0)], accepted=False, reason="criterion")
    run([dict(observed, linear_rhs_norm=0, linear_atol=5e-324,
              linear_residual_limit=5e-324)])
    run([dict(observed, linear_rhs_norm=1e308, linear_rtol=2,
              linear_residual_limit=1e308)], accepted=False, reason="criterion")
    run([dict(observed, linear_residual_limit=0.10000000000000002)])
    run([dict(observed, linear_rhs_norm=0, linear_residual_limit=0.0001)])
    # A failed solve with a valid diagnostic criterion is still readable.
    failed = run([dict(observed, linear_final=10)])
    assert failed["complete"] and failed["loops"][0]["linear_contraction"] == [20.0]
    run([observed, observed])
    # Both ranks can be individually legal while describing different criteria.
    for changed in (dict(observed, linear_rhs_norm=200, linear_residual_limit=0.2),
                    dict(observed, linear_atol=0.01),
                    dict(observed, linear_rtol=0.002, linear_residual_limit=0.2),
                    dict(observed, linear_residual_limit=0.10000000000000002),
                    unavailable):
        run([observed, changed], accepted=False, reason="per-rank linear criterion")
    run([observed, base], schema=3, accepted=False, reason="per-rank loop schema")


def main():
    script = sys.argv[1]
    with tempfile.TemporaryDirectory(prefix="hundun-observe-") as directory:
        root = Path(directory)
        base = dict(step=7, rank=0, attempt=1, corrector=1, refinement=0,
                    kind=0, dropped_loops=0, invoked=1, iterations=2,
                    A_calls=3, M_calls=2, A_ns=30, M_ns=20,
                    baseline_candidates=1, extrapolated_candidates=0,
                    ladder_candidates=0, incomplete_candidates=0,
                    linear_initial=1, linear_final=0.01, globalization_valid=0)
        total = dict(step=7, rank=0, dropped_loops=0, pressure_calls=1,
                     diagonal_calls=0, spatial_calls=0, A_apply_ns=30,
                     M_apply_ns=20)
        (root / "conservation.csv").write_text("step,rank\n7,0\n")
        metadata = ("HUNDUN_V04_THIN_DOMAIN_RUN_V1\nstarting_step 6\n"
                    "requested_steps 1\nexpected_ranks 1\nobservation_schema 3\n"
                    "candidate_identity " + "a" * 64 + "\nend\n")
        (root / "RUN.meta").write_text(metadata)
        identity = hashlib.sha256(metadata.encode()).hexdigest()
        base['source_meta_sha256'] = total['source_meta_sha256'] = identity

        checks = []

        def run(loops, performance, accepted, reason=None):
            write(root / "solver-rank-0.csv", loops)
            write(root / "performance.csv", performance)
            result = subprocess.run([sys.executable, script, str(root)],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, timeout=10)
            assert (result.returncode == 0) == accepted, result.stderr or result.stdout
            if reason:
                assert reason in result.stderr, result.stderr
            checks.append('accepted' if accepted else 'rejected')
            return json.loads(result.stdout) if accepted else None

        old = run([base], [total], True)
        assert old["loops"][0]["scalar_coupling_sweep"] == 1
        assert old["rank_step_mean_ns"]["scalar_remap_ns"] == 0
        check_linear_criterion(script, checks, root / "criterion", base, total)
        loops = [dict(base, scalar_coupling_sweep=i) for i in (1, 2)]
        summed = dict(total, pressure_calls=2, A_apply_ns=60, M_apply_ns=40,
                      scalar_remap_ns=9)
        new = run(loops, [summed], True)
        assert len(new["loops"]) == 2
        assert new["rank_step_mean_ns"]["scalar_remap_ns"] == 9
        run([loops[0], loops[0]], [summed], False)  # Duplicate logical loop.
        run([dict(base, dropped_loops=1)], [total], False)
        run([base], [dict(total, A_apply_ns=29)], False)
        run([base], [dict(total, pressure_calls=0, spatial_calls=1)], False)
        run([base, dict(base, step=8)], [total], False)  # Unaccounted step.
        run([base], [total, total], False)  # Duplicate performance row.
        # Both inputs can be incomplete in exactly the same way. They cannot
        # establish expected rank coverage by comparing with each other.
        metadata = metadata.replace('expected_ranks 1', 'expected_ranks 2')
        (root / 'RUN.meta').write_text(metadata)
        identity = hashlib.sha256(metadata.encode()).hexdigest()
        base['source_meta_sha256'] = total['source_meta_sha256'] = identity
        run([base], [total], False, 'frozen rank set')
        other = dict(base, rank=1)
        write(root / 'solver-rank-1.csv', [other])
        complete = run([base], [total, dict(total, rank=1)], True)
        assert complete['complete'] and complete['expected_ranks'] == [0, 1]
        assert complete['validated_steps'] == 1

        base['scalar_coupling_sweep'] = other['scalar_coupling_sweep'] = 1
        write(root / 'solver-rank-1.csv', [other])
        loops = [dict(base, scalar_coupling_sweep=i) for i in (1, 2)]
        # Keep each rank's own totals internally consistent. The absent second
        # sweep on rank 1 must still fail the independent logical-loop check.
        summed = dict(total, pressure_calls=2, A_apply_ns=60, M_apply_ns=40)
        run(loops, [summed, dict(total, rank=1)], False, 'logical loop')
        write(root / 'solver-rank-1.csv', [other, other])
        run([base], [total, dict(summed, rank=1)], False)
        write(root / 'solver-rank-1.csv', [other])
        for key, value in [('iterations', -1), ('A_ns', -1), ('M_ns', 'nan'),
                           ('linear_final', 'inf'), ('attempt', 0), ('invoked', 2)]:
            run([dict(base, **{key: value})], [total, dict(total, rank=1)], False)
        run([base], [dict(total, final_momentum_ns=-1), dict(total, rank=1)], False)
        run([dict(base, source_meta_sha256='b'*64)], [total, dict(total, rank=1)],
            False, 'source identity mismatch')
        run([base], [dict(total, source_meta_sha256='c'*64), dict(total, rank=1)],
            False, 'source identity mismatch')

        # A fully present-looking last row without its newline is not complete.
        run([base], [total, dict(total, rank=1)], True)
        path = root / 'solver-rank-1.csv'
        path.write_bytes(path.read_bytes().rstrip(b'\r\n'))
        result = subprocess.run([sys.executable, script, str(root)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        assert result.returncode != 0 and 'truncated CSV tail' in result.stderr
        checks.append('truncated tail rejected')
        write(path, [other])

        metadata = metadata.replace('requested_steps 1', 'requested_steps 2')
        (root / 'RUN.meta').write_text(metadata)
        identity = hashlib.sha256(metadata.encode()).hexdigest()
        base[SOURCE] = total[SOURCE] = other[SOURCE] = identity
        write(path, [other])
        run([base], [total, dict(total, rank=1)], False, 'expected step 8')
        partial = subprocess.run([sys.executable, script, str(root), '--allow-partial'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        assert partial.returncode == 0, partial.stderr
        partial = json.loads(partial.stdout)
        assert not partial['complete'] and partial['validated_steps'] == 1 and partial['issues']
        checks.append('partial is not complete')

        final_loop = dict(base, step=8)
        write(path, [other, dict(other, step=8)])
        result = run([base, final_loop], [total, dict(total, rank=1),
            dict(total, step=8), dict(total, step=8, rank=1)], True)
        sidecar = root / 'details.jsonl'
        streamed = subprocess.run([sys.executable, script, str(root), '--details-output', str(sidecar)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        assert streamed.returncode == 0, streamed.stderr
        streamed = json.loads(streamed.stdout)
        assert streamed['complete'] and streamed['loops'] == []
        assert [json.loads(line) for line in sidecar.read_text().splitlines()] == result['loops']
        assert streamed['details']['sha256'] == hashlib.sha256(sidecar.read_bytes()).hexdigest()
        checks.append('streamed details match in-memory result')
        # Two different input sizes, one run each; this is a postprocessor
        # allocation regression, not a repeated CFD timing benchmark.
        peaks = []
        for count in (8, 2048):
            large_meta = metadata.replace('requested_steps 2', 'requested_steps ' + str(count))
            (root / 'RUN.meta').write_text(large_meta)
            identity = hashlib.sha256(large_meta.encode()).hexdigest()
            for rank in (0, 1):
                write(root / ('solver-rank-{}.csv'.format(rank)),
                      [dict(base, step=s, rank=rank, source_meta_sha256=identity)
                       for s in range(7, 7 + count)])
            write(root / 'performance.csv',
                  [dict(total, step=s, rank=rank, source_meta_sha256=identity)
                   for s in range(7, 7 + count) for rank in (0, 1)])
            wrapper = ('import runpy,sys,tracemalloc; tracemalloc.start(); '
                       'path=sys.argv.pop(1); runpy.run_path(path,run_name="__main__"); '
                       'print(tracemalloc.get_traced_memory()[1])')
            measured = subprocess.run([sys.executable, '-c', wrapper, script, str(root),
                '--output', str(root / ('summary-{}.json'.format(count))),
                '--details-output', str(root / ('details-{}.jsonl'.format(count)))],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=20)
            assert measured.returncode == 0, measured.stderr
            peaks.append(int(measured.stdout.strip()))
        assert peaks[1] < peaks[0] + 2 * 1024 * 1024, peaks
        print('streamed CLI traced peak bytes (8/2048 steps): {}'.format(peaks))
        checks.append('bounded streamed Python allocations')
    print("solver observation CLI PASS ({} checks: frozen ranks/range, source identity, corruption, partial, streaming)".format(len(checks)))


if __name__ == "__main__":
    main()
