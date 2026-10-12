"""Mathematical preflight for polynomial evaluation and rejected candidates."""
import importlib.util
import itertools
import json
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('fluid_stencil',Path(__file__).parents[1]/'scripts/audit-m2a-fluid-stencil.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


class FluidStencilTests(unittest.TestCase):
    def test_polynomial_at_origin_and_first_moments(self):
        points=list(itertools.product((-1.5,-.5,.5,1.5),repeat=3))
        weights,defect=audit.weights(points,3)
        polynomial=lambda p:5+2*p[0]-3*p[1]+p[0]*p[2]+4*p[1]**2+7*p[2]**3
        self.assertAlmostEqual(sum(w*polynomial(p) for w,p in zip(weights,points)),5,places=11)
        self.assertLess(defect,1e-10)
        for a in range(3):self.assertAlmostEqual(sum(w*p[a] for w,p in zip(weights,points)),0,places=11)

    def test_three_parallel_layers_cannot_determine_cubic(self):
        points=list(itertools.product((-1.5,-.5,.5,1.5),(-1.5,-.5,.5,1.5),(1.,2.,3.)))
        with self.assertRaisesRegex(ValueError,'Deficient'):audit.weights(points,3)

    def test_known_geometric_candidates_rejected_before_native(self):
        result=audit.run()
        self.assertTrue(all(not candidate['admitted_to_native'] for candidate in result['candidates'][:2]))
        self.assertTrue(result['candidates'][2]['admitted_to_native'])
        self.assertTrue(all(case['min_normalised_cholesky_pivot']>1e-8 for case in result['candidates'][2]['cases']))
        quadratic=result['candidates'][0]['cases'][0]
        self.assertTrue(quadratic['checks']['force_exchange'])
        self.assertFalse(quadratic['checks']['fixed_wall_improvement'])
        self.assertTrue(any('axis_levels' in case.get('error','') for case in result['candidates'][1]['cases']))
        json.dumps(result,allow_nan=False)


if __name__=='__main__':unittest.main()
