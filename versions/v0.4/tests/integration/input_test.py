#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Focused source-layout and provenance checks, independent of MPI/Cantera."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('reference_input', sys.argv.pop(1))
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)
FIXTURE = Path(__file__).parent.parent / 'fixtures' / 'input.d'


class InputTest(unittest.TestCase):
    def setUp(self):
        self.text = FIXTURE.read_text()

    def test_runtime_core_and_provenance(self):
        for profile in ('gtmc', '624cf'):
            result = reader.parse_text(self.text, profile)
            values = result['settings']
            self.assertEqual(values['steps'], 6000)
            self.assertEqual(values['restart_regrid_replace_mesh'], False)
            self.assertEqual(values['statistics'], [False, False, 0])
            self.assertEqual(values['cfl_window'], [.2, .3, .25])
            self.assertEqual(values['reference_properties'], [1.18, .663, 1e-5, 1e5])
            self.assertEqual(values['combustion'], [True, True, 'vode'])
            self.assertEqual(values['pdf_fields'], 4)
            self.assertEqual(values['phase_frequency'], [400., .287275])
            self.assertEqual(values['effective_boundary'], [.28, .005266])
            fields = {f['name']: f for f in result['fields']}
            self.assertEqual(fields['cfl_window']['line'], 37)
            self.assertEqual(fields['reference_properties']['line'], 45)
            for field in fields.values():
                self.assertEqual(field['raw'], self.text.splitlines()[field['line']-1])
            self.assertEqual(result['conversion_status'],
                             'source_catalog_ready_native_case_binding_pending')

    def test_optional_records_and_legacy_layout(self):
        text = '\n'.join(line for line in self.text.splitlines()
                         if 'restart_regrid_replace_mesh' not in line)
        text = text.replace('false false 0 /turbstat', 'false false /turbstat')
        catalog = reader.parse_text(text, 'gtmc')
        runtime = catalog['settings']
        self.assertTrue(runtime['restart_regrid_replace_mesh'])
        self.assertEqual(runtime['statistics'], [False, False, 0])
        self.assertEqual([d['name'] for d in catalog['reader_defaults']],
                         ['restart_regrid_replace_mesh', 'statistics_start_step'])
        legacy = reader.parse_text(text, 'dyn711')['settings']
        self.assertIsNone(legacy['restart_regrid_replace_mesh'])
        self.assertEqual(legacy['statistics'], [False, False])
        self.assertEqual(legacy['pdf_fields'], 4)
        self.assertEqual(legacy['effective_boundary'], [.28, .005266])
        with self.assertRaises(reader.InputError):
            reader.parse_text(self.text, 'dyn711')

    def test_statistics_fallback_and_negative_start(self):
        text = self.text.replace('false false 0 /turbstat', 'false false word /turbstat')
        result = reader.parse_text(text, 'gtmc')
        self.assertEqual(result['settings']['statistics'], [False, False, 0])
        field = next(f for f in result['fields'] if f['name'] == 'statistics')
        self.assertEqual(field['trailing_tokens'], ['word'])
        with self.assertRaisesRegex(reader.InputError, 'nonnegative'):
            reader.parse_text(self.text.replace('false false 0 /turbstat',
                                               'false false -1 /turbstat'), 'gtmc')

    def test_extension_preservation(self):
        result = reader.parse_text(self.text, 'gtmc')
        extensions = result['extensions']
        self.assertEqual(len(extensions), 3)
        self.assertEqual(extensions[0]['values'], ['manual', '2.0'])
        self.assertEqual(extensions[1]['name'], 'time_control')
        self.assertIn('dt_cap = 9.9d-6', extensions[1]['raw'])
        self.assertEqual(extensions[2]['values'], ['name/with!punctuation'])
        with self.assertRaisesRegex(reader.InputError, 'unterminated namelist'):
            reader.parse_text(self.text + '\n&unfinished\n x=1\n', 'gtmc')
        with self.assertRaisesRegex(reader.InputError, 'own record'):
            reader.parse_text(self.text + '\n&inline x=1 /\n', 'gtmc')

    def test_pressure_and_asymmetric_cfl(self):
        result = reader.parse_text(self.text.replace('true                            /compressible',
                                                    'false                           /compressible')
                                  .replace('0.2 0.3 0.25', '0.2 0.4 0.25'), 'gtmc')
        mappings = {m['target']: m for m in result['mappings']}
        self.assertEqual(mappings['flow.thermodynamic_pressure_pa']['value'], 1e5)
        self.assertEqual(mappings['time.convective_cfl_margin']['status'],
                         'asymmetric_window_requires_policy')
        self.assertIsNone(mappings['time.convective_cfl_margin']['value'])

    def test_scalar_syntax(self):
        self.assertEqual(reader.tokens(".true., 2, 'a/b!c', 'it''s' / ignored"),
                         ['.true.', '2', 'a/b!c', "it's"])
        self.assertEqual(reader.real('1.25D-3'), .00125)
        for bad in ('NaN', 'inf', '1D999'):
            with self.assertRaises(reader.InputError): reader.real(bad)
        for bad in (',2', '1,,2', '2*1', "'open"):
            with self.assertRaises(reader.InputError): reader.tokens(bad)
        with self.assertRaises(reader.InputError): reader.logical('...true...')

    def test_rejection_paths(self):
        for text in (self.text[:50], self.text.replace('0.2 0.3 0.25', '0.3 0.2 0.25'),
                     self.text.replace('4 /number of stochastic fields', '3 /number of stochastic fields')):
            with self.assertRaises(reader.InputError): reader.parse_text(text, 'gtmc')
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.d'
            code = Path(directory) / 'input.F90'
            source.write_text(self.text)
            code.write_text('unregistered reader')
            with self.assertRaisesRegex(reader.InputError, 'SHA-256'):
                reader.catalog(source, code)


if __name__ == '__main__':
    unittest.main()
