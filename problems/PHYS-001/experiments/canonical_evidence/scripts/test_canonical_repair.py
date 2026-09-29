"""Bounded regressions for canonical evidence repair; never run the full baseline."""
from __future__ import annotations

import hashlib
import argparse
import contextlib
import io
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import reproduce_and_compare as reproduction
import s2_local_exponent as exponent


class FreshReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.script = self.root / 'canonical_evidence/reference/reference.py'
        self.script.parent.mkdir(parents=True)
        self.script.write_bytes(b'# inert test fixture\n')
        self.output = self.root / 'results/r1/k41_baseline_results.json'
        self.output.parent.mkdir(parents=True)
        self.output.write_text('{"stale": true}', encoding='utf-8')
        self.patchers = [
            mock.patch.object(reproduction, 'REFERENCE_SCRIPT', str(self.script)),
            mock.patch.object(reproduction, 'EXPECTED_REFERENCE_SHA256_LF',
                              hashlib.sha256(self.script.read_bytes()).hexdigest()),
        ]
        for patcher in self.patchers:
            patcher.start()
            self.addCleanup(patcher.stop)

    def test_failed_stale_deletion_never_launches_or_loads_old_json(self):
        with mock.patch.object(reproduction.os, 'remove', side_effect=PermissionError('locked')):
            with mock.patch.object(reproduction.subprocess, 'run') as launch:
                launch.return_value = subprocess.CompletedProcess([], 1, '', '')
                result = reproduction.glm_reference_run()
        launch.assert_not_called()
        self.assertFalse(result['available'])
        self.assertFalse(result.get('results_fresh', False))
        self.assertIsNone(result.get('results'))
        self.assertIn('cannot remove', result['reason'])
        self.assertEqual(json.loads(self.output.read_text()), {'stale': True})

    def test_scientific_fail_with_new_valid_output_is_usable(self):
        def scientific_fail(*args, **kwargs):
            self.assertFalse(self.output.exists())
            self.output.write_text('{"scientific_verdict": "FAIL"}', encoding='utf-8')
            return subprocess.CompletedProcess([], 1, 'expected frozen scientific FAIL', '')
        with mock.patch.object(reproduction.subprocess, 'run', side_effect=scientific_fail):
            result = reproduction.glm_reference_run()
        self.assertTrue(result['results_fresh'])
        self.assertEqual(result['exit'], 1)
        self.assertEqual(result['results']['scientific_verdict'], 'FAIL')

    def test_crash_without_new_output_cannot_reuse_old_json(self):
        with mock.patch.object(reproduction.subprocess, 'run',
                               return_value=subprocess.CompletedProcess([], 1, 'crash', '')):
            result = reproduction.glm_reference_run()
        self.assertIsNone(result['results'])
        self.assertTrue(result['results_stale_or_missing'])
        self.assertFalse(result.get('results_fresh', False))

    def test_malformed_new_output_is_not_fresh_evidence(self):
        def malformed(*args, **kwargs):
            self.output.write_text('partial {', encoding='utf-8')
            return subprocess.CompletedProcess([], 1, '', '')
        with mock.patch.object(reproduction.subprocess, 'run', side_effect=malformed):
            result = reproduction.glm_reference_run()
        self.assertIsNone(result['results'])
        self.assertFalse(result.get('results_fresh', False))

    def test_hash_mismatch_never_launches(self):
        self.script.write_text('# tampered\n', encoding='utf-8')
        with mock.patch.object(reproduction.subprocess, 'run') as launch:
            result = reproduction.glm_reference_run()
        launch.assert_not_called()
        self.assertFalse(result['available'])
        self.assertIn('hash mismatch', result['reason'])


class CutoffInterpretationTests(unittest.TestCase):
    @staticmethod
    def analytic_derivative(r):
        # Independent derivative of the finite sum, using fsum and a stable
        # 2*sin(x/2)**2 denominator rather than subtracting nearby log values.
        frequencies = [2 * math.pi * k / 4096 for k in range(3, 901)]
        numerator = math.fsum(k**-3 * r*w * math.sin(w*r)
                              for k, w in zip(range(3, 901), frequencies))
        denominator = math.fsum(k**-3 * 2*math.sin(w*r/2)**2
                                for k, w in zip(range(3, 901), frequencies))
        return numerator / denominator

    def test_small_lag_near_two_fails_uv_gate(self):
        for r in (1, 2):
            with self.subTest(r=r):
                exact = self.analytic_derivative(r)
                self.assertGreater(exact, 1.9)
                self.assertLess(exponent.K_HI*r, 20)
                self.assertLessEqual(exponent.K_LO*r, 0.1)
                self.assertLess(abs(exponent.local_exponent(r, 3)-exact), 2e-4)

    def test_marginal_window_passes_both_gates_but_is_not_near_two(self):
        for r in (15, 21):
            with self.subTest(r=r):
                exact = self.analytic_derivative(r)
                self.assertGreaterEqual(exponent.K_HI*r, 20)
                self.assertLessEqual(exponent.K_LO*r, 0.1)
                self.assertGreater(exact, 1.70)
                self.assertLess(exact, 1.75)
                self.assertLess(abs(exponent.local_exponent(r, 3)-exact), 2e-4)


class VerdictTests(unittest.TestCase):
    CHECKS = ('agree_k41', 'agree_neg', 'spread_is_zero_k41', 'spread_is_zero_neg')

    def outcome(self, checks, **changes):
        values = dict(gap_k41=2.3e-15, gap_steep=4.6e-15,
                      spread_k41=1.2e-15, spread_steep=3.4e-15,
                      bound=1e-12, spread_bound=1e-12)
        values.update(changes)
        return reproduction.comparison_outcome(checks=checks, **values)

    def test_any_single_failed_check_forbids_success_claim(self):
        for failed in self.CHECKS:
            with self.subTest(failed=failed):
                checks = {name: name != failed for name in self.CHECKS}
                outcome = self.outcome(checks)
                self.assertFalse(outcome['passed'])
                self.assertFalse(outcome['implementation_error_excluded'])
                self.assertEqual(outcome['status'], 'NOT_REPRODUCED')
                self.assertIn(failed, outcome['failed_checks'])
                self.assertIn('Implementation error (B) is not excluded', outcome['verdict'])
                self.assertNotIn('REPRODUCED_WITHIN_FLOAT_BOUND', outcome['verdict'])
                self.assertNotIn('slopes match', outcome['verdict'])

    def test_success_uses_actual_values_not_historical_constants(self):
        outcome = self.outcome({name: True for name in self.CHECKS})
        self.assertTrue(outcome['passed'])
        for value in ('2.3e-15', '4.6e-15', '1.2e-15', '3.4e-15'):
            self.assertIn(value, outcome['verdict'])
        self.assertNotIn('1.1e-14', outcome['verdict'])
        self.assertNotIn('1.3e-14', outcome['verdict'])

    def test_missing_check_and_nonfinite_measurement_fail_closed(self):
        self.assertFalse(self.outcome({})['passed'])
        checks = {name: True for name in self.CHECKS}
        self.assertFalse(self.outcome(checks, gap_k41=None)['passed'])
        self.assertFalse(self.outcome(checks, spread_steep=float('nan'))['passed'])

    def test_main_persists_failure_claim_when_comparison_fails(self):
        reference = {'available': True, 'sha256': 'fixture', 'exit': 1, 'seconds': 0,
                     'results': {'verdict': 'FAIL', 'cases': [
                         {'beta_name': name, 'sf_mean_slope': 1.0,
                          'sf_std': 0.0, 'spec_mean_slope': 1.0}
                         for name in ('K41_positive', 'steep_negative')]}}
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / 'comparison.json'
            with mock.patch.object(reproduction, 'glm_reference_run', return_value=reference), \
                 mock.patch.object(reproduction, 'make_field_own', return_value=[1.0]), \
                 mock.patch.object(reproduction, 's2_own', return_value=1.0), \
                 mock.patch.object(reproduction, 's2_glm_exact', return_value=1.0), \
                 mock.patch.object(reproduction, 'ols', return_value=0.0), \
                 mock.patch.object(reproduction, 'slope_over', return_value=0.0), \
                 contextlib.redirect_stdout(io.StringIO()) as stdout:
                status = reproduction.main(['--json', str(destination)])
            comparison = json.loads(destination.read_text(encoding='utf-8'))['comparison']
        self.assertEqual(status, 1)
        self.assertEqual(comparison['status'], 'NOT_REPRODUCED')
        self.assertFalse(comparison['implementation_error_excluded'])
        self.assertNotIn('REPRODUCED_WITHIN_FLOAT_BOUND:', stdout.getvalue())
        self.assertIn('Implementation error (B) is not excluded', stdout.getvalue())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path)
    parser.add_argument('--log', type=Path)
    args = parser.parse_args()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__]))
    log = stream.getvalue()
    print(log, end='')
    if args.log:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        args.log.write_text(log, encoding='utf-8')
    if args.json:
        rows = []
        for r in (1, 2, 15, 21):
            analytic = CutoffInterpretationTests.analytic_derivative(r)
            finite_difference = exponent.local_exponent(r, 3)
            rows.append({'r': r, 'analytic_derivative': analytic,
                         'finite_difference': finite_difference,
                         'absolute_gap': abs(analytic-finite_difference),
                         'k_lo_r': exponent.K_LO*r, 'k_hi_r': exponent.K_HI*r,
                         'ir_gate': exponent.K_LO*r <= .1, 'uv_gate': exponent.K_HI*r >= 20})
        report = {'scope': 'software repair and fixed-point finite-sum regression only',
                  'historical_scientific_evidence_regenerated': False,
                  'cutoff_rows': rows,
                  'regression_tests': {'count': result.testsRun,
                                       'failures': len(result.failures),
                                       'errors': len(result.errors),
                                       'passed': result.wasSuccessful()},
                  'environment': {'python': sys.version, 'platform': platform.platform()},
                  'original_scientific_states': {'E3': 'FAIL', 'B4': 'NOT_MET_AS_FROZEN',
                                                 'all_criteria_passed': False, 'PHYS-001': 'OPEN'}}
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    raise SystemExit(0 if result.wasSuccessful() else 1)
