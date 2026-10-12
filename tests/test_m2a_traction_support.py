import importlib.util
from pathlib import Path
import sys
import unittest

root = Path(__file__).parents[1]
sys.path.insert(0, str(root/'scripts'))
spec = importlib.util.spec_from_file_location('traction_support', root/'scripts/audit-m2a-traction-support.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class TractionSupport(unittest.TestCase):
    def test_weighted_overlap_distinguishes_zero_weight_corner(self):
        self.assertEqual(audit.overlap({(0,): 0., (1,): .25, (2,): .75}, {(0,), (1,)}), .25)

    def test_known_geometry_and_source_free_far_probe(self):
        result = audit.fixture(64, 0.)
        self.assertEqual(result['support_cells'], 776)
        self.assertGreater(result['exterior_support_cells'], 0)
        self.assertEqual(result['traction_probes_distance_h']['4']['nodes_with_support_overlap'], 0)


if __name__ == '__main__':
    unittest.main()
