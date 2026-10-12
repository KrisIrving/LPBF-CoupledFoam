"""Additional static FV ledger identities; never replace original physical gates."""
import csv
import math


def evaluate(case, meta):
    controls = meta['control_volume_audit']
    with (case/'control-volume-history.csv').open(newline='') as stream:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    with (case/'mechanical-history.csv').open(newline='') as stream:
        times = [float(row['time']) for row in csv.DictReader(stream)]
    if not rows or any(not math.isfinite(v) for row in rows for v in row.values()):
        raise ValueError('Missing/nonfinite control volume ledger')
    if len(times) != meta['expected_rows'] or len(rows) != len(times)*len(controls['padding_h']):
        raise ValueError('Incomplete control volume trajectory')
    for window, time in enumerate(times):
        pair = rows[window*len(controls['padding_h']):(window+1)*len(controls['padding_h'])]
        if [row['padding_h'] for row in pair] != controls['padding_h'] or any(abs(row['time']-time) > 1e-12 for row in pair):
            raise ValueError('Control volume time/padding mismatch or duplicate')
    limits = {'matrix_match_impulse': controls['matrix_match_impulse_N_s'],
              'face_match_impulse': controls['face_match_impulse_N_s'],
              'angular_face_match_impulse': controls['angular_match_impulse_N_m_s'],
              'source_force_error': controls['source_force_N'],
              'source_torque_error': controls['source_torque_N_m']}
    maxima = {key: max(row[key] for row in rows) for key in limits}
    checks = {key: all(0 <= row[key] <= limit for row in rows) for key, limit in limits.items()}
    checks['positive_cell_count'] = all(row['cells'] > 0 and row['cells'].is_integer() for row in rows)
    linear = angular = 0.
    for row in rows:
        linear = max(linear, math.sqrt(sum((row['inventory_'+a]-row['source_'+a]-row['face_net_'+a]-row['residual_'+a])**2 for a in 'xyz'))*meta['delta_t'])
        angular = max(angular, math.sqrt(sum((row['inventory_torque_'+a]-row['source_torque_'+a]-row['face_torque_'+a]-row['arm_correction_'+a]-row['residual_torque_'+a])**2 for a in 'xyz'))*meta['delta_t'])
    checks['linear_identity'] = linear <= controls['face_match_impulse_N_s']
    checks['angular_identity_with_arm_correction'] = angular <= controls['angular_match_impulse_N_m_s']
    maxima.update(linear_identity_impulse=linear, angular_identity_impulse=angular,
                  cell_impulse_L1=max(row['cell_impulse_L1'] for row in rows))
    # The old full-equation gate already limits L1; this diagnostic does not
    # use a cancelling control-volume sum to supersede it.
    return {'passed': all(checks.values()), 'checks': checks, 'metrics': maxima,
            'rows': len(rows), 'scope': __doc__}
