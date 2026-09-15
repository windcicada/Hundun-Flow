#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exercise public algorithm selection, evidence and cross-rank Restart."""
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def main():
    executable, launcher, validator = map(lambda p: str(Path(p).resolve()), sys.argv[1:4])
    assert sys.argv[4:] in ([], ['--cantera']), sys.argv[4:]
    cantera = sys.argv[4:] == ['--cantera']

    def command(args, success=True):
        result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                universal_newlines=True, timeout=45)
        if (result.returncode == 0) != success:
            raise RuntimeError('{}\n{}'.format(' '.join(args), result.stdout[-8000:]))
        return result.stdout

    def native(ranks, *args, **kwargs):
        return command([launcher, '-n', str(ranks), executable] + [str(a) for a in args], **kwargs)

    def check(case, time, coupling):
        text = native(1, 'check', case, '--dry-plan')
        assert 'time_scheme='+time in text and 'coupling='+coupling in text, text
        return re.search(r'VALID case=(\d+) product=(\d+)', text).groups()

    def evidence(output, time, coupling, previous=0, manifest=None):
        path = output / 'evidence.jsonl'
        args = [sys.executable, validator, 'runtime', str(path)]
        if manifest:
            args += ['--run-start-manifest', str(manifest)]
        command(args)
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        assert len(rows) == 1 and rows[0]['step'] == previous + 1
        row = rows[0]
        assert row['algorithm'] == {'time_scheme': time, 'coupling': coupling}
        assert row['run_start']['previous_step'] == previous
        if coupling == 'outer_corrected':
            assert row['cold']['outer_iterations'] >= 1 and row['momentum_predictor_passes'] == 0
        else:
            assert row['momentum_predictor_passes'] == (1 if coupling == 'PISO' else 2)

    with tempfile.TemporaryDirectory(prefix='hf-alg-') as directory:
        root = Path(directory)
        case = root / 'c'
        native(1, 'init-case', '--output', case)
        source = json.loads((case / 'case.json').read_text())
        assert source['time']['scheme'] == 'cn_be'
        assert source['solver']['coupling'] == 'outer_corrected'

        def variant(name, time, coupling):
            target = root / name
            target.mkdir()
            model = json.loads(json.dumps(source))
            model['time']['scheme'] = time
            model['solver']['coupling'] = coupling
            (target / 'case.json').write_text(json.dumps(model))
            (target / 'thermophysics.d').write_bytes((case / 'thermophysics.d').read_bytes())
            return target

        old = variant('old', 'cn_be', 'PISO')
        assert check(case, 'cn_be', 'outer_corrected') == check(old, 'cn_be', 'outer_corrected')
        common = ['--steps', '1', '--max-dt', '1e-5', '--output-interval', '0',
                  '--diagnostics-interval', '0']
        seed = root / 'seed'
        # A historical input writes the source; the canonical input reads it
        # on four ranks with the original method history intact.
        native(1, 'run', old, '--output', seed, *common, '--restart-interval', '1')
        evidence(seed, 'cn_be', 'outer_corrected')
        manifests = list((seed / 'Restart').glob('generation-*/manifest.bin'))
        assert len(manifests) == 1
        resumed = root / 'next'
        native(4, 'run', case, '--output', resumed, *common,
               '--restart', seed / 'Restart', '--restart-interval', '0')
        evidence(resumed, 'cn_be', 'outer_corrected', 1, manifests[0])
        for name, coupling in [('p', 'PISO'), ('s', 'SIMPLE')]:
            model = variant(name, 'backward_euler', coupling)
            check(model, 'backward_euler', coupling)
            output = root / (name+'out')
            native(2, 'run', model, '--output', output, *common, '--restart-interval', '0')
            evidence(output, 'backward_euler', coupling)
        invalid = variant('bad', 'backward_euler', 'outer_corrected')
        native(1, 'check', invalid, success=False)

        # The native analytic provider carries distinct formation enthalpies:
        # this exercises heat release and EOS closure as well as source wiring.
        fixture = Path(__file__).resolve().parents[1] / 'fixtures' / 'reacting-isomer'
        reactive = root / 'r'
        reactive.mkdir()
        model = json.loads((fixture / 'case.json').read_text())
        model['time']['scheme'] = 'cn_be'
        model['solver']['coupling'] = 'outer_corrected'
        model['solver']['cold_stopping'] = dict(reference_time=1., momentum=1e-7,
                                               enthalpy=1e-10, species=1e-10)
        model['flow']['pressure_reference'] = 'boundary_absolute'
        for side in ('x_min', 'x_max'):
            boundary = model['boundaries'][side]
            outlet = side == 'x_max'
            boundary['flow_kind'] = 'pressure_outlet' if outlet else 'no_slip_wall'
            boundary['thermal_kind'] = 'none' if outlet else 'adiabatic_wall'
            boundary['scalars'] = [dict(stable_name='A', kind='zero_gradient', value=0,
                backflow_kind='dirichlet' if outlet else 'zero_gradient',
                backflow_value=.25 if outlet else 0)]
            if outlet:
                boundary.update(pressure=101325, allow_backflow=True, backflow_temperature=300)
        (reactive / 'case.json').write_text(json.dumps(model))
        (reactive / 'thermophysics.d').write_bytes((fixture / 'thermophysics.d').read_bytes())
        check(reactive, 'cn_be', 'outer_corrected')
        reaction_options = ['--steps', '1', '--output-interval', '0', '--diagnostics-interval', '0']
        seed = root / 'r1'
        text = native(1, 'run', reactive, '--output', seed, *reaction_options,
                      '--initial-state', '101325,300,0,0,0,0.25', '--restart-interval', '1')
        assert 'source=transport_reactor' in text
        assert 'chemistry=transport_then_interval' in native(1, 'check', reactive)
        evidence(seed, 'cn_be', 'outer_corrected')
        manifests = list((seed / 'Restart').glob('generation-*/manifest.bin'))
        assert len(manifests) == 1
        resumed = root / 'r4'
        native(4, 'run', reactive, '--output', resumed, *reaction_options,
               '--restart', seed / 'Restart', '--restart-interval', '0')
        evidence(resumed, 'cn_be', 'outer_corrected', 1, manifests[0])

        # Explicit SGS coefficient participates in the compiled model and Restart.
        smag = root / 'sm'
        smag.mkdir()
        smodel = json.loads((reactive / 'case.json').read_text())
        smodel['turbulence'] = dict(model='smagorinsky', coefficient=.12)
        (smag / 'case.json').write_text(json.dumps(smodel))
        (smag / 'thermophysics.d').write_bytes((fixture / 'thermophysics.d').read_bytes())
        description = native(1, 'check', smag)
        assert 'sgs=smagorinsky' in description and 'smagorinsky_coefficient=0.12' in description
        seed, resumed = root / 'sm1', root / 'sm4'
        native(1, 'run', smag, '--output', seed, *reaction_options,
               '--initial-state', '101325,300,0,0,0,0.25', '--restart-interval', '1')
        evidence(seed, 'cn_be', 'outer_corrected')
        manifests = list((seed / 'Restart').glob('generation-*/manifest.bin'))
        assert len(manifests) == 1
        native(4, 'run', smag, '--output', resumed, *reaction_options,
               '--restart', seed / 'Restart', '--restart-interval', '0')
        evidence(resumed, 'cn_be', 'outer_corrected', 1, manifests[0])
        smodel['turbulence']['coefficient'] = .18
        (smag / 'case.json').write_text(json.dumps(smodel))
        native(1, 'run', smag, '--output', root / 'sm_bad', *reaction_options,
               '--restart', seed / 'Restart', '--restart-interval', '0', success=False)

        # Repeat the native boundary/Restart path with PaSR and active WALE.
        pasr = root / 'rp'
        pasr.mkdir()
        model['reaction']['model'] = 'pasr_algebraic_v1'
        model['reaction']['mixing'] = dict(c_z=.001, turbulent_schmidt=.7)
        model['turbulence']['model'] = 'wale'
        (pasr / 'case.json').write_text(json.dumps(model))
        (pasr / 'thermophysics.d').write_bytes((fixture / 'thermophysics.d').read_bytes())
        check(pasr, 'cn_be', 'outer_corrected')
        seed = root / 'p1'
        text = native(1, 'run', pasr, '--output', seed, *reaction_options,
                      '--initial-state', '101325,300,0,0,0,0.25', '--restart-interval', '1')
        assert 'source=frozen_pasr' in text
        evidence(seed, 'cn_be', 'outer_corrected')
        manifests = list((seed / 'Restart').glob('generation-*/manifest.bin'))
        assert len(manifests) == 1
        resumed = root / 'p4'
        native(4, 'run', pasr, '--output', resumed, *reaction_options,
               '--restart', seed / 'Restart', '--restart-interval', '0')
        evidence(resumed, 'cn_be', 'outer_corrected', 1, manifests[0])
        if cantera:
            fixture = Path(__file__).resolve().parents[1] / 'fixtures' / 'reacting-cantera-isomer'
            provider = json.loads((fixture / 'case.json').read_text())['reaction']
            for name, mode in [('cm', 'finite_rate_mean'), ('cp', 'pasr_algebraic_v1')]:
                case = root / name
                case.mkdir()
                model['reaction'] = dict(provider, model=mode)
                if mode == 'pasr_algebraic_v1':
                    model['reaction']['mixing'] = dict(c_z=.001, turbulent_schmidt=.7)
                (case / 'case.json').write_text(json.dumps(model))
                for asset in ('thermophysics.d', 'mechanism.yaml'):
                    (case / asset).write_bytes((fixture / asset).read_bytes())
                check(case, 'cn_be', 'outer_corrected')
                description = native(1, 'check', case)
                response_rtol = float(re.search(r'source_response_rtol=([^ ]+)', description).group(1))
                response_atol = float(re.search(r'source_response_atol=([^ ]+)', description).group(1))
                if mode == 'finite_rate_mean':
                    assert 'energy=conservative_total' in description
                    assert 0 < response_rtol <= 1.000001e-11
                    assert 0 < response_atol <= response_rtol
                else:
                    assert response_rtol == response_atol == 0

                seed, resumed = root / (name+'1'), root / (name+'4')
                native(1, 'run', case, '--output', seed, *reaction_options,
                       '--initial-state', '101325,300,0,0,0,0.25', '--restart-interval', '1')
                evidence(seed, 'cn_be', 'outer_corrected')
                manifests = list((seed / 'Restart').glob('generation-*/manifest.bin'))
                assert len(manifests) == 1
                native(4, 'run', case, '--output', resumed, *reaction_options,
                       '--restart', seed / 'Restart', '--restart-interval', '0')
                evidence(resumed, 'cn_be', 'outer_corrected', 1, manifests[0])
    print('algorithm CLI: cold/mean/PaSR CN/BE 1->4 Restart, BE PISO/SIMPLE; Cantera={}'.format(cantera))


if __name__ == '__main__':
    main()
