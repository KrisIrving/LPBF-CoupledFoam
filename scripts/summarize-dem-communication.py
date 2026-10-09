"""Independent checks of user-side embedded DEM feedback, not CFD validation."""
import argparse
import csv
import json
import math
from pathlib import Path


def check_case(path, ranks):
    try:
        if (path / 'result.txt').read_text().strip() != 'exit_status=0':
            raise ValueError('Case execution failed')
        if f'DEMO_PASS ranks={ranks} windows=20' not in (path / 'log.demo').read_text():
            raise ValueError('Missing success marker')
        with (path / 'communication-history.csv').open(newline='') as f:
            rows = list(csv.DictReader(f))
        if len(rows) != 40:
            raise ValueError('Expected 20 windows and 2 particles')
        expected_mass = 2500 * 4 / 3 * math.pi * 0.005**3
        previous = {1: ([0.45, 0.3, 0.5], [0.4, 0., 0.]),
                    2: ([0.55, 0.7, 0.5], [-0.4, 0., 0.])}
        owners = {}
        states = {}
        migrations = 0
        max_x = max_v = max_impulse = max_clock = 0.
        seen = set()
        for row in rows:
            values = {k: float(v) for k, v in row.items()}
            if not all(math.isfinite(v) for v in values.values()):
                raise ValueError('Nonfinite diagnostic')
            w, pid, owner = (int(values[k]) for k in ('window', 'id', 'owner'))
            if any(values[k] != n for k, n in zip(('window', 'id', 'owner'), (w, pid, owner))):
                raise ValueError('Noninteger window, ID or owner')
            if pid not in (1, 2) or w not in range(1, 21) or owner not in range(ranks):
                raise ValueError('Invalid window, ID or rank')
            if (w, pid) in seen or (w > 1 and (w-1, pid) not in seen):
                raise ValueError('Duplicate or out-of-order particle history')
            seen.add((w, pid))
            mass = values['mass']
            if abs(mass / expected_mass - 1) > 1e-12 or abs(values['radius'] - 0.005) > 1e-14:
                raise ValueError('Particle mass or radius changed')
            x0, v0 = previous[pid]
            force = [values[k] for k in ('fx', 'fy', 'fz')]
            target = 0.6 if pid == 1 else -0.6
            expected_force = [mass * 2 * (target-v0[0]), 0., 0.]
            if any(abs(a-b) > 1e-14 for a, b in zip(force, expected_force)):
                raise ValueError('Host feedback force does not follow measured velocity')
            x = [values[k] for k in ('x', 'y', 'z')]
            v = [values[k] for k in ('vx', 'vy', 'vz')]
            for k in range(3):
                a = expected_force[k] / mass
                max_x = max(max_x, abs(x[k] - x0[k] - v0[k]*0.01 - 0.5*a*0.01**2))
                max_v = max(max_v, abs(v[k] - v0[k] - a*0.01))
                max_impulse = max(max_impulse, abs(mass*(v[k]-v0[k])-expected_force[k]*0.01))
            max_clock = max(max_clock, abs(values['time_of'] - w*0.01),
                            abs(values['time_dem'] - w*0.01))
            if pid in owners and owners[pid] != owner:
                migrations += 1
            owners[pid] = owner
            previous[pid] = x, v
            states[(w, pid)] = x + v
        if max_x > 1e-10 or max_v > 1e-10 or max_impulse > 1e-12 or max_clock > 1e-12:
            raise ValueError('Independent motion, impulse or clock check failed')
        if ranks == 2 and migrations < 2:
            raise ValueError('Both particles must demonstrate DEM-rank migration')
        return {'passed': True, 'windows': 20, 'migrations_observed': migrations,
                'max_position_error_m': max_x, 'max_velocity_error_m_s': max_v,
                'max_impulse_residual_kg_m_s': max_impulse, 'max_clock_error_s': max_clock}, states
    except (OSError, ValueError, KeyError, TypeError) as e:
        return {'passed': False, 'reason': str(e)}, {}


def summarize(root):
    serial, s = check_case(root / 'serial', 1)
    mpi, m = check_case(root / 'mpi', 2)
    comparison = {'passed': False, 'reason': 'Both trajectories must pass first'}
    if serial['passed'] and mpi['passed']:
        difference = max(abs(a-b) for k in s for a, b in zip(s[k], m[k]))
        comparison = {'passed': difference <= 1e-10, 'max_state_absolute_difference': difference,
                      'position_tolerance_m': 1e-10, 'velocity_tolerance_m_s': 1e-10}
    return {'scope': 'OpenFOAM host and embedded LIGGGHTS communication only; no fluid solve, '
                     'resolved particle forces, heat, restart or scalable distributed CFD mapping.',
            'serial': serial, 'mpi': mpi, 'comparison': comparison,
            'passed': serial['passed'] and mpi['passed'] and comparison['passed']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('report', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.report)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result['passed'] else 1)
