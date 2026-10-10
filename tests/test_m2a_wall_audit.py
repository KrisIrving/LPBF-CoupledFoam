import importlib.util
from pathlib import Path
import sys
import unittest

scripts=Path(__file__).parents[1]/'scripts';sys.path.insert(0,str(scripts))
spec=importlib.util.spec_from_file_location('wall_audit',scripts/'audit-m2a-wall-operators.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


class WallAuditTests(unittest.TestCase):
    def test_fine_analytic_stress_approaches_independent_exact_load(self):
        fixed=audit.fixture('fixed',256,0.,'linear')
        self.assertLess(abs(fixed['pressure_drag_ratio']-1/3),.01)
        self.assertLess(abs(fixed['viscous_drag_ratio']-2/3),.01)
        rotate=audit.fixture('rotate',256,0.,'linear')
        self.assertLess(abs(rotate['stress_load_ratio']-1),.01)

    def test_ghost_order_does_not_change_fluid_only_stress_oracle(self):
        a=audit.fixture('fixed',36,.25,'linear')
        b=audit.fixture('fixed',36,.25,'quadratic')
        self.assertEqual(a['pressure_force'],b['pressure_force'])
        self.assertEqual(a['viscous_force'],b['viscous_force'])
        self.assertNotEqual(a['interpolated_image_wall_slip_relative'],b['interpolated_image_wall_slip_relative'])


if __name__=='__main__':unittest.main()
