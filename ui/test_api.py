"""No production solver is launched by these API tests."""
import json
import os
import tempfile
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch, Mock

from fastapi.testclient import TestClient

from ui.backend import create_app
from ui.jobs import Registry, ApiError, process
from ui.data import ident, rows


class ApiTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.case = self.root/'cases'/'c'
        self.case.mkdir(parents=True)
        (self.case/'case.json').write_text(json.dumps({'mesh': {'exact_cells': [4, 3, 2]}, 'time': {'scheme': 'cn_be'}}))
        self.run = self.case/'run'
        self.run.mkdir()
        (self.run/'status.json').write_text('{"phase":"solving","step":8,"time":0.5}')
        (self.run/'monitor.jsonl').write_text('{"step":8,"time":0.5,"payload":{"dt":0.1,"seconds":2,"continuity":1e-9}}\n{"step":9')
        self.registry = Registry([self.root/'cases'], self.root/'state', self.root/'missing', self.root/'missing-mpi')
        self.client = TestClient(create_app(self.registry), base_url='http://localhost')
        self.key = ident(self.run)

    def tearDown(self):
        self.temp.cleanup()

    def test_real_rows_and_stale_status(self):
        result = self.client.get('/api/runs').json()['runs'][0]
        self.assertEqual(result['step'], 8)
        self.assertEqual(result['status'], 'unverified')
        self.assertIsNone(result['alive'])
        self.assertEqual(result['process_status'], 'unknown')
        payload = self.client.get('/api/runs/'+self.key).json()
        self.assertEqual(len(payload['history']), 1)
        self.assertEqual(payload['history'][0]['seconds'], 2)
        self.assertEqual(payload['case']['mesh']['exact_cells'], [4, 3, 2])

    def test_root_and_origin_guards(self):
        response = self.client.post('/api/runs/import', json={'path': '/etc'})
        self.assertEqual(response.status_code, 400)
        response = self.client.post('/api/runs/import', json={'path': str(self.run)}, headers={'Origin': 'http://evil.example'})
        self.assertEqual(response.status_code, 403)
        self.assertEqual(self.client.get('/api/health', headers={'Host': 'evil.example'}).status_code, 403)

    def test_evidence_fills_unsampled_monitor_and_preserves_units(self):
        row = {'step': 9, 'time': .6, 'max_rank_step_ns': 3000000000,
            'max_rank_rss_bytes': 2048, 'terminal_physical_audit': {
                'continuity_residual': 2e-9, 'energy_residual': 3e-8,
                'committed_convective_cfl': {'dt': .1, 'out_max': .2}},
            'cold': {'outer_iterations': 3, 'pressure_iterations': 10}}
        (self.run/'evidence.jsonl').write_text(json.dumps(row)+'\n')
        payload = self.client.get('/api/runs/'+self.key).json()
        self.assertEqual([r['step'] for r in payload['history']], [8, 9])
        self.assertEqual(payload['history'][-1]['seconds'], 3)
        self.assertEqual(payload['history'][-1]['cfl'], .2)
        self.assertEqual(payload['resources']['max_rank_rss_bytes'], 2048)

    def test_missing_measurement_remains_null(self):
        (self.run/'monitor.jsonl').unlink()
        (self.run/'evidence.jsonl').write_text('{"step":1,"time":0.1}\n')
        payload = self.client.get('/api/runs/'+self.key).json()
        self.assertIsNone(payload['history'][0]['seconds'])
        self.assertIsNone(payload['history'][0]['continuity'])

    def test_bad_run_does_not_break_list(self):
        (self.run/'status.json').write_text('[]')
        result = self.client.get('/api/runs')
        self.assertEqual(result.status_code, 200)
        self.assertEqual(result.json()['runs'][0]['status'], 'unavailable')

    def test_published_fields_only_run_is_read_only_with_unknown_clock(self):
        path = self.root/'cases'/'fields'
        visit = path/'Visit'
        visit.mkdir(parents=True)
        (visit/'step-00000001.rank0.vtk').write_text('# vtk DataFile Version 3.0\n')
        (visit/'step-00000001.visit').write_text('!NBLOCKS 1\nstep-00000001.rank0.vtk\n')
        self.registry.refresh(force=True)
        item = self.registry.summary(ident(path))
        self.assertIsNone(item['step'])
        self.assertIsNone(item['time'])
        self.assertTrue(item['has_fields'])
        self.assertFalse(any(item['capabilities'].values()))

    def test_log_symlinks_and_special_files_are_not_read(self):
        secret = self.root/'secret'
        secret.write_text('PRIVATE CONTENT')
        (self.run/'stdout.log').symlink_to(secret)
        os.mkfifo(str(self.run/'stderr.log'))
        payload = self.client.get('/api/runs/'+self.key).json()
        self.assertEqual(payload['logs'], [])

    def test_finished_control_rejected(self):
        response = self.client.post('/api/runs/'+self.key+'/control', json={'action': 'pause', 'request_id': 'stop-123456'})
        self.assertEqual(response.status_code, 409)
        self.assertFalse((self.run/'stop').exists())

    def test_external_output_process_is_live_but_read_only(self):
        info = {'pid': 123, 'starttime': 'same', 'argv': ['/external/hundun', 'run', str(self.case), '--output', str(self.run)], 'cwd': self.root, 'mpi': 128}
        self.registry._processes = [info]
        self.registry._process_time = time.monotonic()
        with patch('ui.jobs.process', return_value=info):
            item = self.registry.summary(self.key)
            self.assertTrue(item['alive'])
            self.assertEqual(item['process_status'], 'external')
            self.assertEqual(item['mpi'], 128)
            self.assertEqual(item['status'], 'solving')
            self.assertFalse(any(item['capabilities'].values()))
            with self.assertRaises(ApiError):
                self.registry.control(self.key, 'pause', 'external-stop-01')
            with self.assertRaises(ApiError):
                self.registry.launch(self.key, 1, 1, 'external-resume-01', resume=True)
        self.assertFalse((self.run/'stop').exists())

    def test_exit_evidence_and_recorded_mpi(self):
        (self.run/'status.json').write_text('{"phase":"completed","step":8}')
        (self.run/'evidence.jsonl').write_text('{"step":8,"mpi_size":4}\n')
        item = self.registry.summary(self.key)
        self.assertFalse(item['alive'])
        self.assertEqual(item['process_status'], 'exited')
        self.assertEqual(item['status'], 'completed')
        self.assertEqual(item['mpi'], 4)

    def test_output_is_not_checkpoint_and_idempotent(self):
        with patch.object(self.registry, 'live', return_value={'pid': 7}):
            first = self.registry.control(self.key, 'output', 'out-123456')
            self.assertTrue((self.run/'output').is_file())
            self.assertFalse((self.run/'Restart').exists())
            (self.run/'output').unlink()
            self.assertEqual(first, self.registry.control(self.key, 'output', 'out-123456'))
            self.assertFalse((self.run/'output').exists())
            with self.assertRaises(ApiError):
                self.registry.control(self.key, 'pause', 'out-123456')

    def test_control_restart_receipt_semantics(self):
        with patch.object(self.registry, 'live', return_value={'pid': 7}):
            result = self.registry.control(self.key, 'pause', 'pause-123456')
        self.assertTrue((self.run/'stop').is_file())
        self.assertEqual(result['control']['published'], 'Restart/current')
        self.assertTrue(result['pending'])

    def test_launch_idempotency_persisted_before_dispatch(self):
        calls = []
        def operation(db):
            calls.append(1)
            return {'ok': True}
        self.registry.once('request-123456', ['launch'], operation)
        other = Registry(self.registry.roots, self.registry.state, self.registry.program, self.registry.mpi)
        other.once('request-123456', ['launch'], operation)
        self.assertEqual(len(calls), 1)

    def test_parallel_duplicate_request_dispatches_once(self):
        calls = []
        def operation(db):
            calls.append(1)
            time.sleep(.03)
            return {'ok': True}
        with ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(lambda _: self.registry.once('parallel-123456', ['start'], operation), range(2)))
        self.assertEqual(len(calls), 1)
        self.assertEqual(results[0], results[1])

    def test_pid_reuse_and_path_identity(self):
        self.registry.program = Path('/safe/run')
        info = {'argv': ['/safe/bin/hundun', 'run', str(self.case), '--output', str(self.run)], 'cwd': self.root}
        self.assertTrue(self.registry.native_process(info, self.run))
        info['argv'][0] = '/unregistered/bin/hundun'
        self.assertFalse(self.registry.native_process(info, self.run))
        info['argv'][0] = '/safe/bin/hundun'
        info['argv'][-1] = '/another/run'
        self.assertFalse(self.registry.native_process(info, self.run))
        db = self.registry.database()
        db['jobs'][self.key] = {'path': str(self.run), 'pid': 123, 'starttime': 'old'}
        self.registry.save(db)
        with patch('ui.jobs.process', return_value={'starttime': 'new'}):
            self.assertIsNone(self.registry.live(self.key))

    def test_numbers_and_paths_never_become_commands(self):
        response = self.client.post('/api/cases/'+ident(self.case)+'/start', json={'steps': '1;touch /tmp/no', 'ranks': 1, 'request_id': 'test-123456'})
        self.assertEqual(response.status_code, 422)
        response = self.client.post('/api/cases/'+ident(self.case)+'/start', json={'steps': 1, 'ranks': 0, 'request_id': 'test-123457'})
        self.assertEqual(response.status_code, 400)

    def test_start_uses_new_directory_and_recovers_registration(self):
        self.registry.program.write_text('#!/bin/sh\nexit 0\n')
        self.registry.program.chmod(0o700)
        self.registry.mpi.write_text('#!/bin/sh\nexit 0\n')
        self.registry.mpi.chmod(0o700)
        with patch.object(self.registry, 'command', return_value='VALID'), patch('ui.jobs.subprocess.Popen', return_value=Mock(pid=os.getpid())) as launch:
            result = self.registry.launch(ident(self.case), 2, 1, 'launch-123456')
            repeated = self.registry.launch(ident(self.case), 2, 1, 'launch-123456')
            with self.assertRaises(ApiError) as failure:
                self.registry.launch(ident(self.case), 2, 1, 'launch-654321')
            self.assertEqual(failure.exception.code, 'case_running')
        self.assertEqual(result, repeated)
        self.assertEqual(launch.call_count, 1)
        self.assertNotEqual(result['run']['path'], str(self.run))
        argv = launch.call_args[0][0]
        self.assertIn('--output', argv)
        self.assertNotIn('shell', launch.call_args[1])
        recovered = Registry(self.registry.roots, self.registry.state, self.registry.program, self.registry.mpi)
        self.assertIn(result['run']['id'], recovered.runs)
        self.assertEqual(recovered.database()['jobs'][result['run']['id']]['pid'], os.getpid())

    def test_import_persists_and_symlink_escape_rejected(self):
        result = self.registry.import_run(str(self.run))
        self.registry.refresh()
        self.assertIn(result['id'], self.registry.runs)
        outside = self.root/'outside'
        outside.mkdir()
        (outside/'status.json').write_text('{}')
        (self.case/'escape').symlink_to(outside)
        with self.assertRaises(ApiError):
            self.registry.import_run(str(self.case/'escape'))

    def test_registered_run_cannot_be_retargeted_with_symlink(self):
        moved = self.case/'moved'
        self.run.rename(moved)
        self.run.symlink_to(moved)
        with self.assertRaises(ApiError) as failure:
            self.registry.resolve_run(self.key)
        self.assertEqual(failure.exception.code, 'path_changed')

    def test_resume_rejects_mismatched_case_identity_before_launch(self):
        (self.run/'status.json').write_text('{"phase":"completed","step":8}')
        self.registry.program.write_text('#!/bin/sh\nexit 0\n')
        self.registry.program.chmod(0o700)
        self.registry.mpi.write_text('#!/bin/sh\nexit 0\n')
        self.registry.mpi.chmod(0o700)
        (self.run/'Restart').mkdir()
        (self.run/'Restart/current').write_text('generation')
        (self.run/'evidence.jsonl').write_text('{"step":8,"case":13}\n')
        with patch.object(self.registry, 'command', return_value='VALID case=12 product=99'), patch('ui.jobs.subprocess.Popen') as launch:
            with self.assertRaises(ApiError) as failure:
                self.registry.launch(self.key, 2, 1, 'resume-123456', resume=True)
        self.assertEqual(failure.exception.code, 'case_identity_mismatch')
        self.assertFalse(launch.called)


if __name__ == '__main__':
    unittest.main()
