"""Offline shell scheduling checks; no OpenFOAM, DEM or WSL invocation."""
import os
from pathlib import Path
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
BASH = os.environ.get('LPBF_TEST_BASH') or shutil.which('bash')


@unittest.skipUnless(BASH, 'Bash unavailable')
class CaseQueueTests(unittest.TestCase):
    def run_queue(self, limit, failing='none'):
        script = '''
set -eo pipefail
source scripts/case-queue.sh
worker() {
    printf 'start %s\\n' "$1"
    sleep .1
    printf 'end %s\\n' "$1"
    [ "$1" != "$failing" ]
}
if lpbf_run_case_queue "$limit" worker a b c d e; then rc=0; else rc=$?; fi
printf 'queue_exit %s\\n' "$rc"
'''
        result = subprocess.run([BASH, '-c', script], cwd=ROOT, text=True,
                                capture_output=True, timeout=15,
                                env={**os.environ, 'limit': str(limit), 'failing': failing})
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.splitlines()

    def check_events(self, lines, limit):
        active = set()
        starts, ends, peak = set(), set(), 0
        for line in lines:
            event, value = line.split()
            if event == 'start':
                self.assertNotIn(value, starts)
                starts.add(value)
                active.add(value)
                peak = max(peak, len(active))
                self.assertLessEqual(len(active), limit)
            elif event == 'end':
                self.assertIn(value, active)
                active.remove(value)
                ends.add(value)
            else:
                self.assertFalse(active, 'Queue returned before workers completed')
        self.assertEqual(starts, set('abcde'))
        self.assertEqual(ends, starts)
        return peak

    def test_bounded_concurrency(self):
        lines = self.run_queue(2)
        self.assertEqual(self.check_events(lines, 2), 2)
        self.assertEqual(lines[-1], 'queue_exit 0')

    def test_failure_preserves_all_remaining_workers(self):
        lines = self.run_queue(2, 'b')
        self.check_events(lines, 2)
        self.assertEqual(lines[-1], 'queue_exit 1')

    def test_serial_and_invalid_limit(self):
        lines = self.run_queue(1)
        self.assertEqual(self.check_events(lines, 1), 1)
        self.assertEqual(self.run_queue(0), ['queue_exit 2'])


if __name__ == '__main__':
    unittest.main()
