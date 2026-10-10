"""Independent fixture checks against known legacy geometry and analytic volume."""
import importlib.util
from pathlib import Path
import sys
import unittest

SCRIPTS = Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
spec = importlib.util.spec_from_file_location('geometry_audit',SCRIPTS/'m2a-geometry-audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class GeometryMatrixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = audit.run()

    def test_six_fixtures_and_analytic_closure(self):
        self.assertEqual(len(self.result['cases']),6)
        self.assertTrue(self.result['passed'])
        for case in self.result['cases'].values():
            self.assertLess(case['metrics']['volume_relative_error'],1e-6)
            self.assertEqual(case['metrics']['shared_face_relative_difference'],0)

    def test_independent_legacy_metrics_and_phase_sensitivity(self):
        cases = self.result['cases']
        self.assertAlmostEqual(cases['n36-phase0.0']['metrics']['old_q4_support_volume_ratio'],
                               1.4147106052612897,places=12)
        self.assertAlmostEqual(cases['n48-phase0.0']['metrics']['old_q4_volume_relative_error'],
                               .005753551511382016,places=12)
        self.assertGreater(cases['n36-phase0.25']['metrics']['old_q4_support_volume_ratio'],
                           cases['n36-phase0.0']['metrics']['old_q4_support_volume_ratio'])


if __name__ == '__main__':
    unittest.main()
