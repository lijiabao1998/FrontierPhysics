"""Replay the original seven-test harness against the immutable pre-fix runner.

The only expected failing test is stale-output deletion. No baseline is executed;
the harness mocks reference subprocesses. Git blobs are exported to a temporary
directory inside this checkout's ignored scratch directory; no branch is changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[5]
PREFIX = 'problems/PHYS-001/experiments/canonical_evidence/scripts/'
SUBJECT_SHA = 'e3c93deb54384c93e039cd3330602d0a0d08c604'
HARNESS_SHA = 'e68111cb11f2b2e777e1195e525904345f321472'


def main():
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path, required=True)
    parser.add_argument('--log', type=Path, required=True)
    args = parser.parse_args()
    scratch = ROOT / 'scratch'
    scratch.mkdir(exist_ok=True)
    exports = [('reproduce_and_compare.py', SUBJECT_SHA),
               ('s2_local_exponent.py', SUBJECT_SHA),
               ('test_canonical_repair.py', HARNESS_SHA)]
    provenance = []
    with tempfile.TemporaryDirectory(prefix='prefixed-control-', dir=scratch) as directory:
        target = Path(directory).resolve()
        if target.parent != scratch.resolve():
            raise RuntimeError('temporary export escaped the checkout scratch directory')
        for name, commit in exports:
            blob = subprocess.check_output(['git', '-C', str(ROOT), 'show', commit+':'+PREFIX+name])
            (target / name).write_bytes(blob)
            provenance.append({'path': PREFIX+name, 'commit': commit,
                               'sha256_git_blob': hashlib.sha256(blob).hexdigest()})
        command = [sys.executable, 'test_canonical_repair.py']
        process = subprocess.run(command, cwd=target, capture_output=True,
                                 text=True, encoding='utf-8', errors='replace',
                                 env={**os.environ, 'PYTHONIOENCODING': 'utf-8'})
    combined = process.stdout + process.stderr
    expected = (process.returncode == 1 and 'Ran 7 tests' in combined
                and 'FAILED (failures=1)' in combined
                and 'FAIL: test_failed_stale_deletion_never_launches_or_loads_old_json' in combined
                and 'Expected \'run\' to not have been called. Called 1 times.' in combined)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.log.write_text(combined, encoding='utf-8')
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps({'subject_sha': SUBJECT_SHA, 'harness_sha': HARNESS_SHA,
                                    'exports': provenance, 'command': command,
                                    'working_directory': 'isolated export under repository scratch/',
                                    'observed_exit': process.returncode,
                                    'expected_negative_control_observed': expected}, indent=2)+'\n',
                         encoding='utf-8')
    print(combined, end='')
    print('Expected pre-fix negative control observed:', expected)
    return 0 if expected else 1


if __name__ == '__main__':
    raise SystemExit(main())
