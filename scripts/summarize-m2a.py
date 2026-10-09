"""Independent mechanical ledgers and predeclared M2A-01 gates; no CFD execution."""
import argparse
import csv
import json
import math
from pathlib import Path

CASES = ('fixed-coarse', 'fixed-fine', 'rotate-coarse', 'rotate-fine',
         'translate', 'free', 'free-restart', 'free-mpi2')


def vector(row, names):
    return [row[k] for k in names]


def norm(values):
    return math.sqrt(sum(x*x for x in values))


def columns(prefix):
    return [prefix+c for c in 'xyz']


def failure_diagnostics(case):
    """Retain incomplete trajectory evidence without evaluating it as a pass."""
    try:
        with (case/'mechanical-history.csv').open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        result = {'recorded_windows': len(rows), 'completion': False}
        if rows:
            for key in ('time', 'vx', 'fx', 'tz', 'slip_rms', 'momentum_residual'):
                value = float(rows[-1][key])
                result['last_'+key] = value if math.isfinite(value) else None
        return result
    except (OSError, ValueError, KeyError, TypeError):
        return {'completion': False, 'trajectory_diagnostics': 'unavailable'}


def evaluate(case):
    meta = json.loads((case/'M2A_META.json').read_text())
    limits = meta['thresholds']
    if (case/'result.txt').read_text().strip() != 'exit_status=0':
        raise ValueError('Case execution failed; inspect retained logs')
    logs = ('log.first', 'log.restart') if meta['selection'] == 'free-restart' else ('log.solver',)
    for log in logs:
        if 'M2A_EXECUTION_COMPLETE' not in (case/log).read_text(errors='replace'):
            raise ValueError(f'Missing completion marker: {log}')
    with (case/'mechanical-history.csv').open(newline='') as stream:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream)]
    if len(rows) != meta['expected_rows']:
        raise ValueError(f'Expected {meta["expected_rows"]} windows, got {len(rows)}')
    if any(not math.isfinite(v) for r in rows for v in r.values()):
        raise ValueError('Nonfinite CSV value')
    dt, radius = meta['delta_t'], meta['radius']
    mass = 4*math.pi*radius**3*meta['particle_density']/3
    inertia = 0.4*mass*radius**2
    scale = (meta['reference_flow_speed'] if meta['mode'] == 'fixed' else
             meta['reference_rotation_speed']*radius if meta['mode'] == 'rotate' else
             norm(meta['initial_velocity']))
    checks, metrics = {}, {}
    maxima = {k: 0. for k in ('momentum', 'angular', 'force_decomposition',
                             'torque_decomposition', 'fluid_continuity')}
    previous_v, previous_w = meta['initial_velocity'], meta['initial_omega']
    previous_fluid = None
    ranks = 2 if meta['selection'] == 'free-mpi2' else 1
    for index, row in enumerate(rows, 1):
        if abs(row['time']-index*dt) > limits['clock_s']:
            raise ValueError('Missing, duplicate or out-of-order coupling window')
        if row['id'] != meta['id'] or row['owner'] != int(row['owner']) or not 0 <= row['owner'] < ranks:
            raise ValueError('Invalid stable ID/owner')
        velocity, omega = vector(row, columns('v')), vector(row, columns('w'))
        force, torque = vector(row, columns('f')), vector(row, columns('t'))
        new = vector(row, columns('fluid_p'))
        old = vector(row, columns('old_fluid_p'))
        boundary = vector(row, columns('boundary_f'))
        constraint = vector(row, columns('constraint_f'))
        correction = vector(row, columns('inertia_f'))
        ct, it = vector(row, columns('constraint_t')), vector(row, columns('inertia_t'))
        support = [0., 0., 0.] if meta['mode'] == 'free' else [-f for f in force]
        residual = norm([new[k]-old[k]+mass*(velocity[k]-previous_v[k])
                         -dt*(boundary[k]+support[k]) for k in range(3)])
        maxima['momentum'] = max(maxima['momentum'], residual)
        if meta['mode'] == 'free':
            maxima['angular'] = max(maxima['angular'], norm([
                inertia*(omega[k]-previous_w[k])-dt*torque[k] for k in range(3)]))
        else:
            if norm([velocity[k]-meta['initial_velocity'][k] for k in range(3)]) > 1e-12:
                raise ValueError('Prescribed translation changed')
            if norm([omega[k]-meta['initial_omega'][k] for k in range(3)]) > 1e-12:
                raise ValueError('Prescribed rotation changed')
        maxima['force_decomposition'] = max(maxima['force_decomposition'],
            norm([force[k]+constraint[k]-correction[k] for k in range(3)]))
        maxima['torque_decomposition'] = max(maxima['torque_decomposition'],
            norm([torque[k]+ct[k]-it[k] for k in range(3)]))
        if previous_fluid is not None:
            maxima['fluid_continuity'] = max(maxima['fluid_continuity'],
                                            norm([old[k]-previous_fluid[k] for k in range(3)]))
        previous_v, previous_w, previous_fluid = velocity, omega, new
    mapping = {'momentum': 'momentum_residual_kg_m_s',
               'angular': 'angular_impulse_residual_kg_m2_s',
               'force_decomposition': 'force_decomposition_N',
               'torque_decomposition': 'torque_decomposition_N_m',
               'fluid_continuity': 'momentum_residual_kg_m_s'}
    for key, threshold in mapping.items():
        metrics[key] = maxima[key]
        checks[key] = maxima[key] <= limits[threshold]
    for key, column, threshold, divisor in (
        ('volume', 'volume_error', 'volume_relative', 1.),
        ('slip', 'slip_rms', 'slip_relative', scale),
        ('divergence', 'div_max', 'divergence_per_s', 1.),
        ('clock', 'clock_error', 'clock_s', 1.)):
        metrics[key] = max(abs(r[column])/divisor for r in rows)
        checks[key] = metrics[key] <= limits[threshold]
    if meta['mode'] in ('fixed', 'rotate'):
        exact = (6*math.pi*meta['rho']*meta['nu']*radius*meta['reference_flow_speed']
                 if meta['mode'] == 'fixed' else
                 8*math.pi*meta['rho']*meta['nu']*radius**3*meta['reference_rotation_speed'])
        load = sum(r['fx'] if meta['mode'] == 'fixed' else -r['tz'] for r in rows[-10:])/10
        metrics['stokes_relative_error'] = abs(load/exact-1)
        checks['stokes_reference'] = metrics['stokes_relative_error'] <= limits[
            'stokes_force_relative' if meta['mode'] == 'fixed' else 'stokes_torque_relative']
    if meta['mode'] == 'free':
        checks['actual_response'] = 0 < rows[-1]['vx'] < meta['initial_velocity'][0]-1e-8
    if meta['mode'] == 'translate':
        checks['motion_across_plane'] = meta['initial_centre'][0] < 0 < rows[-1]['x']
    if meta['selection'] == 'free-mpi2':
        checks['multi_rank_coverage'] = all(r['covered_ranks'] == 2 for r in rows)
    if meta['selection'] == 'free-restart':
        checks['common_checkpoint'] = all((case/'0.002'/name).is_file()
                                          and (case/'0.002'/name).stat().st_size > 0
                                          for name in ('couplingState', 'dem.restart'))
    return {'passed': all(checks.values()), 'checks': checks, 'metrics': metrics}, rows, meta


def compare(reference, candidate, limits):
    groups = {'position': (list('xyz'), 'comparison_position_m'),
              'velocity': (columns('v'), 'comparison_velocity_m_s'),
              'rotation': (columns('w'), 'comparison_rotation_per_s'),
              'force': (columns('f'), 'comparison_force_N'),
              'torque': (columns('t'), 'comparison_torque_N_m'),
              'fluid_momentum': (columns('fluid_p'), 'comparison_fluid_momentum_kg_m_s')}
    if len(reference) != len(candidate):
        raise ValueError('Comparison history lengths differ')
    metrics = {group: max(abs(a[col]-b[col]) for a,b in zip(reference,candidate) for col in cols)
               for group, (cols, _) in groups.items()}
    checks = {group: metrics[group] <= limits[threshold] for group, (_,threshold) in groups.items()}
    return {'passed': all(checks.values()), 'checks': checks, 'metrics': metrics}


def summarize(root):
    result = {'package': 'M2A-01', 'cases': {}, 'comparisons': {},
              'physical_validation': 'Pending: experimental mechanical code verification only; '
                                     'not full LPBF or general resolved CFD-DEM validation.'}
    histories, metadata = {}, {}
    for name in CASES:
        try:
            entry, histories[name], metadata[name] = evaluate(root/name)
            result['cases'][name] = entry
        except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError) as error:
            result['cases'][name] = {'passed': False, 'error': str(error),
                                     'diagnostics': failure_diagnostics(root/name)}
    for name in ('free-restart', 'free-mpi2'):
        if 'free' in histories and name in histories:
            result['comparisons'][name] = compare(histories['free'], histories[name], metadata['free']['thresholds'])
        else:
            result['comparisons'][name] = {'passed': False, 'status': 'not_evaluated',
                                           'error': 'Missing completed valid histories'}
    for family in ('fixed', 'rotate'):
        coarse, fine = result['cases'][family+'-coarse'], result['cases'][family+'-fine']
        if 'stokes_relative_error' in coarse.get('metrics', {}) and 'stokes_relative_error' in fine.get('metrics', {}):
            delta = fine['metrics']['stokes_relative_error']-coarse['metrics']['stokes_relative_error']
            result['comparisons'][family+'-refinement'] = {
                'passed': delta <= metadata[family+'-fine']['thresholds']['refinement_error_increase'],
                'error_increase': delta, 'interpretation': 'Two resolutions; no asymptotic convergence claim'}
        else:
            result['comparisons'][family+'-refinement'] = {'passed': False, 'status': 'not_evaluated',
                                                         'error': 'Missing reference metrics'}
    result['passed'] = all(v['passed'] for section in ('cases', 'comparisons') for v in result[section].values())
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('report', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    summary = summarize(args.report)
    args.output.write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    for name, entry in summary['cases'].items():
        print(name, 'PASS' if entry['passed'] else 'FAIL', entry.get('error', ''))
    for name, entry in summary['comparisons'].items():
        print(name, 'PASS' if entry['passed'] else 'FAIL')
    raise SystemExit(0 if summary['passed'] else 1)
