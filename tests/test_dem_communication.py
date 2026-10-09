import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('dem_summary', Path(__file__).resolve().parents[1]
                                           / 'scripts/summarize-dem-communication.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(root, ranks):
    root.mkdir()
    (root / 'result.txt').write_text('exit_status=0\n')
    (root / 'log.demo').write_text(f'DEMO_PASS ranks={ranks} windows=20\n')
    columns = ('window,id,owner,time_of,time_dem,x,y,z,vx,vy,vz,fx,fy,fz,mass,radius,'
               'x_error,v_error,impulse_error,clock_error').split(',')
    m = 2500 * 4 / 3 * module.math.pi * 0.005**3
    x, v = [0.45, 0.55], [0.4, -0.4]
    with (root / 'communication-history.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        for w in range(1, 21):
            for i in range(2):
                a = 2 * ((0.6 if i == 0 else -0.6)-v[i])
                x[i] += v[i]*0.01 + 0.5*a*0.01**2
                v[i] += a*0.01
                owner = int(x[i] >= 0.5) if ranks == 2 else 0
                writer.writerow([w, i+1, owner, w*0.01, w*0.01, x[i], 0.3+0.4*i,
                                 0.5, v[i], 0, 0, m*a, 0, 0, m, 0.005, 0, 0, 0, 0])


class CommunicationContract(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        fixture(self.root / 'serial', 1)
        fixture(self.root / 'mpi', 2)

    def tearDown(self):
        self.tmp.cleanup()

    def modify(self, fn):
        p = self.root / 'mpi/communication-history.csv'
        with p.open(newline='') as f:
            reader = csv.DictReader(f)
            keys, rows = reader.fieldnames, list(reader)
        fn(rows)
        with p.open('w', newline='') as f:
            writer = csv.DictWriter(f, keys)
            writer.writeheader()
            writer.writerows(rows)

    def test_complete_exchange(self):
        result = module.summarize(self.root)
        self.assertTrue(result['passed'])
        self.assertEqual(result['mpi']['migrations_observed'], 2)

    def test_reported_zero_errors_do_not_hide_bad_motion(self):
        self.modify(lambda rows: rows[12].update(vx='0.7'))
        self.assertFalse(module.summarize(self.root)['passed'])

    def test_clock_drift_rejected(self):
        self.modify(lambda rows: rows[-1].update(time_dem='0.199'))
        self.assertFalse(module.summarize(self.root)['passed'])

    def test_force_transfer_rejected(self):
        self.modify(lambda rows: rows[0].update(fx='0'))
        self.assertFalse(module.summarize(self.root)['passed'])

    def test_missing_particle_rejected(self):
        self.modify(lambda rows: rows.pop())
        self.assertFalse(module.summarize(self.root)['passed'])

    def test_stationary_ownership_rejected(self):
        self.modify(lambda rows: [row.update(owner='0') for row in rows])
        self.assertFalse(module.summarize(self.root)['passed'])

    def test_missing_execution_result_rejected(self):
        (self.root / 'mpi/result.txt').unlink()
        self.assertFalse(module.summarize(self.root)['passed'])


if __name__ == '__main__':
    unittest.main()
