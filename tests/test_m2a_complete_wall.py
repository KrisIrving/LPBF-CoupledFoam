"""Rank-selection and paired wall sampling contracts; no CFD runtime."""
import importlib.util
from pathlib import Path
import unittest
import sys

root=Path(__file__).parents[1]
sys.path.insert(0,str(root/'scripts'))
spec=importlib.util.spec_from_file_location('complete_wall',root/'scripts/audit-m2a-complete-wall.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


class CompleteWallTests(unittest.TestCase):
    def test_duplicate_rows_removed_without_span_loss(self):
        selected,pivots,error=audit.select_rows([{0:1.},{1:1.},{0:.5,1:.5},{0:1.}])
        self.assertEqual(len(selected),2)
        self.assertGreater(min(pivots),1e-8)
        self.assertLess(error,1e-12)

    def test_near_dependent_direction_is_reported_not_hidden(self):
        selected,_,error=audit.select_rows([{0:1.},{0:1.,1:1e-5}])
        self.assertEqual(len(selected),1)
        self.assertGreater(error,9e-6)

    def test_cell_sampling_and_transpose_preserve_rigid_work_and_moments(self):
        point=[.007,-.003,.005];origin=[-.058]*3;h=.004
        row=audit.sampling_row(point,origin,h)
        self.assertAlmostEqual(sum(row.values()),1.)
        for a in range(3):
            self.assertAlmostEqual(sum(w*(origin[a]+k[a]*h) for k,w in row.items()),point[a])
        force=.17;field={k:k[0]+2*k[1]-k[2] for k in row}
        grid_work=sum(force*w*field[k] for k,w in row.items())
        self.assertAlmostEqual(grid_work,force*sum(w*field[k] for k,w in row.items()))

    def test_coarse_offset_cannot_claim_complete_wall_coverage(self):
        candidate=audit.fixture(36,.25)
        self.assertFalse(candidate['admitted_to_native'])
        self.assertGreater(candidate['max_unresolved_row_relative_norm'],1e-5)


if __name__=='__main__':unittest.main()
