"""Manufactured ledger gates, including noncentral face-force torque."""
import csv
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

root = Path(__file__).parents[1]
sys.path.insert(0, str(root/'scripts'))
import control_volume_gates as gates
spec = importlib.util.spec_from_file_location('cv_prepare', root/'scripts/prepare-m2a-joint.py')
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class ControlVolumeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.case = Path(self.temp.name)
        self.meta = {'expected_rows': 1, 'delta_t': .1, 'control_volume_audit': {
            'padding_h': [6, 8], 'matrix_match_impulse_N_s': 1e-12,
            'face_match_impulse_N_s': 1e-12, 'angular_match_impulse_N_m_s': 1e-14,
            'source_force_N': 1e-12, 'source_torque_N_m': 1e-14}}
        # A force along y has face moment 3z but cell-centre moment 2z.
        # The arm correction -1z is required, even with exact linear balance.
        self.rows = []
        for padding in (6, 8):
            row = {'time': .1, 'padding_h': padding, 'cells': 2., 'cell_impulse_L1': 0.}
            for key in ('matrix_match_impulse','face_match_impulse','angular_face_match_impulse','source_force_error','source_torque_error'):
                row[key] = 0.
            for name in ('inventory','source','face_net','residual','inventory_torque','source_torque','face_torque','arm_correction','residual_torque'):
                for axis in 'xyz': row[name+'_'+axis] = 0.
            row.update(inventory_y=1., face_net_y=1., inventory_torque_z=2., face_torque_z=3., arm_correction_z=-1.)
            self.rows.append(row)
        (self.case/'mechanical-history.csv').write_text('time\n0.1\n')

    def tearDown(self):
        self.temp.cleanup()

    def evaluate(self):
        with (self.case/'control-volume-history.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(self.rows[0]))
            writer.writeheader(); writer.writerows(self.rows)
        return gates.evaluate(self.case, self.meta)

    def test_noncentral_face_requires_explicit_arm_correction(self):
        self.assertTrue(self.evaluate()['passed'])
        self.rows[0]['arm_correction_z'] = 0.
        self.assertFalse(self.evaluate()['checks']['angular_identity_with_arm_correction'])

    def test_missing_or_duplicate_window_is_not_accepted(self):
        self.rows[1]['padding_h'] = 6.
        with self.assertRaises(ValueError): self.evaluate()
        self.rows.pop()
        with self.assertRaises(ValueError): self.evaluate()

    def test_source_coverage_and_matrix_mismatch_fail(self):
        self.rows[0]['source_force_error'] = 1e-9
        self.rows[1]['matrix_match_impulse'] = 1e-10
        entry = self.evaluate()
        self.assertFalse(entry['checks']['source_force_error'])
        self.assertFalse(entry['checks']['matrix_match_impulse'])

    def test_symmetric_traction_remains_a_separate_diagnostic(self):
        self.meta['control_volume_audit']['symmetric_stress_audit'] = True
        for row in self.rows:
            for name in ('symmetric_extra_force', 'symmetric_extra_torque',
                         'symmetric_face_net', 'symmetric_face_torque'):
                for axis in 'xyz': row[name+'_'+axis] = 0.
            row.update(symmetric_face_net_y=1., symmetric_face_torque_z=3.4,
                       symmetric_extra_torque_z=.4)
        self.assertTrue(self.evaluate()['passed'])
        self.rows[1]['symmetric_extra_torque_z'] = 0.
        self.assertFalse(self.evaluate()['checks']['symmetric_stress_components'])

    def test_opt_in_preserves_original_physical_controls(self):
        meta = prepare.prepare(self.case, 'fixed-coarse', control_volume_audit=True)
        self.assertEqual(meta['control_volume_audit']['padding_h'], [6, 8])
        self.assertEqual(meta['joint_controls']['krylov_budget'], 64)
        self.assertIn('controlVolumeAudit true;', (self.case/'constant/mechanicalProperties').read_text())
        with self.assertRaises(ValueError):
            prepare.prepare(self.case, 'fixed-coarse', spatial_only=True, control_volume_audit=True)

    def test_preflight_stops_after_one_window(self):
        prepare.prepare(self.case, 'fixed-coarse', control_volume_audit=True)
        prepare.one_window(self.case)
        import json
        meta = json.loads((self.case/'M2A_META.json').read_text())
        self.assertEqual(meta['expected_rows'], 1)
        self.assertTrue(meta['preflight_only'])
        self.assertIn('endTime 0.0001;', (self.case/'system/controlDict').read_text())


if __name__ == '__main__': unittest.main()
