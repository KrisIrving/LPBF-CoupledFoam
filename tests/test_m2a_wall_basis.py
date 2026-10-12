"""Independent modal force/torque/work and complete sampling-space oracle."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np

root=Path(__file__).parents[1];sys.path.insert(0,str(root/'scripts'))
spec=importlib.util.spec_from_file_location('modal_oracle',root/'scripts/audit-m2a-wall-basis.py')
oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)


class WallBasisTests(unittest.TestCase):
    def test_all_geometries_have_complete_space_and_exchange(self):
        for n in (36,48,64):
            for phase in (0.,.25):
                with self.subTest(n=n,phase=phase):
                    result=oracle.fixture(n,phase)
                    self.assertTrue(result['admitted_geometry'])
                    self.assertGreater(result['smallest_retained_singular'],1e-7)
                    self.assertLess(result['largest_discarded_singular'],1e-12)
        self.assertEqual(oracle.fixture(36,.25)['modes'],244)

    def test_zero_sum_mode_has_torque_but_no_point_force(self):
        basis=np.array([[1.,1.],[-1.,1.]])/np.sqrt(2)
        position=np.array([[-1.,0.,0.],[1.,0.,0.]])
        multiplier=np.array([[0.,0.,0.],[0.,1.,0.]])
        source=basis.T@multiplier
        sums=basis.sum(axis=1);first=basis@position
        force=(sums[:,None]*multiplier).sum(axis=0)
        torque=np.cross(first,multiplier).sum(axis=0)
        np.testing.assert_allclose(force,[0,0,0],atol=1e-14)
        np.testing.assert_allclose(torque,[0,0,np.sqrt(2)],atol=1e-14)
        np.testing.assert_allclose(torque,np.cross(position,source).sum(axis=0),atol=1e-14)
        # A rigid velocity must transform with moments, including translation.
        v=np.array([.2,.3,.4]);omega=np.array([.1,.2,.3])
        actual=basis@(v+np.cross(omega,position))
        np.testing.assert_allclose(actual,sums[:,None]*v+np.cross(omega,first),atol=1e-14)


if __name__=='__main__':unittest.main()
