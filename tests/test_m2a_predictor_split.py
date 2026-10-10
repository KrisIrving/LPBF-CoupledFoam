"""Structural counterexample, not a CFD reproduction or accuracy oracle.

A is SPD, strictly diagonally dominant (I plus a weighted graph Laplacian).
D u = u0-u1 = 0, J u = u2 = 0, G=D^T, S=J^T.
Both algorithms enforce D/J exactly at every iteration, but replacing the
H-only update by a full predictor while keeping a diagonal wall response can
make the momentum/pressure/load iteration diverge.
"""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('split_dense_oracle',
    Path(__file__).parents[1]/'scripts/constraint_projection_reference.py')
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)
A = [[50., -18., -24., -7.], [-18., 45., -7., -19.],
     [-24., -7., 44., -12.], [-7., -19., -12., 39.]]


def step(state, full):
    u, p, force = state[:4], state[4], state[5]
    source = [-p, p, force, 0.]
    if full:
        u = oracle.solve(A, source)
    h = [( (force if i == 2 else 0.)
          - sum(A[i][j]*u[j] for j in range(4) if i != j))/A[i][i]
         for i in range(4)]
    p = (h[0]-h[1])/(1/A[0][0]+1/A[1][1])
    u = [h[0]-p/A[0][0], h[1]+p/A[1][1], h[2], h[3]]
    # Diagonal pressure-projected wall response: D R S is zero here.
    force -= u[2]*A[2][2]
    u[2] = 0.
    return u+[p, force]


class PredictorSplitTests(unittest.TestCase):
    def trajectory(self, full, count):
        state = [1., 2., 3., 4., 1., 1.]
        for _ in range(count):
            state = step(state, full)
            scale = max(1., max(abs(x) for x in state))
            self.assertLessEqual(abs(state[0]-state[1]), scale*1e-14)
            self.assertEqual(state[2], 0.)
        return max(abs(x) for x in state)

    def test_h_only_converges_with_diagonal_wall_response(self):
        self.assertLess(self.trajectory(False, 100), 1e-12)

    def test_full_predictor_can_diverge_despite_exact_constraints(self):
        self.assertGreater(self.trajectory(True, 12), 1e10)


if __name__ == '__main__':
    unittest.main()
