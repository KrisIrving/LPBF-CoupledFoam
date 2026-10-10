"""Offline analysis/generation tests, not substitutes for native CFD verification."""
import csv
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT/'scripts'/filename)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


prepare = module('prepare_m2a', 'prepare-m2a.py')
audit = module('audit_m2a', 'summarize-m2a.py')


class ContractTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def fixture(self, selection='free'):
        case = self.root/selection
        meta = prepare.prepare(case, selection)
        (case/'result.txt').write_text('exit_status=0\n')
        logs = ['log.first', 'log.restart'] if selection == 'free-restart' else ['log.solver']
        for log in logs:
            (case/log).write_text('M2A_EXECUTION_COMPLETE\n')
        if selection == 'free-restart':
            (case/'0.002').mkdir()
            for name in ['couplingState', 'dem.restart']:
                (case/'0.002'/name).write_text('synthetic test checkpoint')
        rows = []
        mass = 4*math.pi*meta['radius']**3*meta['particle_density']/3
        speed, dt, force = meta['initial_velocity'][0], meta['delta_t'], -0.0001
        fluid, x = 0., meta['initial_centre'][0]
        for i in range(40):
            row = {key: 0. for key in ['time', 'id', 'owner', 'x', 'y', 'z', 'vx', 'vy', 'vz',
                   'wx', 'wy', 'wz', 'fx', 'fy', 'fz', 'tx', 'ty', 'tz', 'volume_error',
                   'slip_rms', 'div_max', 'momentum_residual', 'angular_impulse_residual',
                   'clock_error', 'covered_ranks']}
            row.update(pressure_correctors=8, momentum_equation_impulse_L1=0., support_volume_ratio=1.2)
            for prefix in ['fluid_p', 'old_fluid_p', 'boundary_f', 'constraint_f', 'inertia_f',
                           'constraint_t', 'inertia_t']:
                row.update({prefix+k: 0. for k in 'xyz'})
            row.update(time=(i+1)*dt, id=7, old_fluid_px=fluid, fx=force, constraint_fx=-force,
                       covered_ranks=2 if selection == 'free-mpi2' else 1)
            speed += force*dt/mass
            x += speed*dt
            fluid -= force*dt
            row.update(vx=speed, x=x, fluid_px=fluid)
            rows.append(row)
        self.save(case, rows)
        return case, rows

    @staticmethod
    def save(case, rows):
        with (case/'mechanical-history.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    def test_consistent_sparse_id_ledger(self):
        case, _ = self.fixture()
        result, _, _ = audit.evaluate(case)
        self.assertTrue(result['passed'], result)

    def test_false_native_zero_does_not_hide_bad_fluid_balance(self):
        case, rows = self.fixture()
        rows[8]['fluid_px'] += 1e-6
        self.save(case, rows)
        result, _, _ = audit.evaluate(case)
        self.assertFalse(result['checks']['momentum'])
        self.assertFalse(result['checks']['fluid_continuity'])

    def test_boundary_and_force_are_independently_checked(self):
        case, rows = self.fixture()
        rows[12]['boundary_fx'] = 0.01
        rows[13]['constraint_fx'] = 0.
        self.save(case, rows)
        result, _, _ = audit.evaluate(case)
        self.assertFalse(result['checks']['momentum'])
        self.assertFalse(result['checks']['force_decomposition'])

    def test_duplicate_window_and_nonfinite_rejected(self):
        case, rows = self.fixture()
        rows[4]['time'] = rows[3]['time']
        self.save(case, rows)
        with self.assertRaisesRegex(ValueError, 'window'):
            audit.evaluate(case)
        rows[4]['time'] = 0.0005
        rows[9]['vz'] = float('nan')
        self.save(case, rows)
        with self.assertRaisesRegex(ValueError, 'Nonfinite'):
            audit.evaluate(case)

    def test_restart_checkpoint_and_both_runs_required(self):
        case, _ = self.fixture('free-restart')
        self.assertTrue(audit.evaluate(case)[0]['passed'])
        (case/'0.002/dem.restart').unlink()
        self.assertFalse(audit.evaluate(case)[0]['checks']['common_checkpoint'])
        (case/'log.restart').write_text('failed before completion')
        with self.assertRaisesRegex(ValueError, 'marker'):
            audit.evaluate(case)

    def test_parallel_requires_actual_cfd_coverage(self):
        case, rows = self.fixture('free-mpi2')
        self.assertTrue(audit.evaluate(case)[0]['passed'])
        rows[10]['covered_ranks'] = 1
        self.save(case, rows)
        self.assertFalse(audit.evaluate(case)[0]['checks']['multi_rank_coverage'])

    def test_complete_matrix_and_comparison_difference(self):
        for selection in prepare.CASES:
            meta = prepare.prepare(self.root/selection, selection)
            self.assertLess(0.12/meta['mesh_n'], 2*meta['radius'])
            self.assertEqual(meta['delta_t']/meta['dem_delta_t'], 10)
        _, rows = self.fixture()
        other = [dict(row) for row in rows]
        other[5]['wx'] = 0.01
        limits = prepare.prepare(self.root/'check', 'free')['thresholds']
        self.assertFalse(audit.compare(rows, other, limits)['checks']['rotation'])
        result = audit.summarize(self.root)
        self.assertFalse(result['passed'])
        self.assertEqual(len(result['cases']), 8)

    def test_failed_run_keeps_partial_evidence_without_passing(self):
        case, rows = self.fixture()
        self.save(case, rows[:1])
        (case/'result.txt').write_text('exit_status=1\n')
        summary = audit.summarize(self.root)
        entry = summary['cases']['free']
        self.assertFalse(entry['passed'])
        self.assertEqual(entry['diagnostics']['recorded_windows'], 1)
        self.assertEqual(entry['diagnostics']['last_time'], 0.0001)
        self.assertEqual(summary['comparisons']['free-restart']['status'], 'not_evaluated')

    def test_global_balance_does_not_hide_unconverged_local_equations(self):
        case, rows = self.fixture()
        rows[0]['momentum_equation_impulse_L1'] = 1e-8
        self.save(case, rows)
        result, _, _ = audit.evaluate(case)
        self.assertTrue(result['checks']['momentum'])
        self.assertFalse(result['checks']['equation_convergence'])


if __name__ == '__main__':
    unittest.main()
