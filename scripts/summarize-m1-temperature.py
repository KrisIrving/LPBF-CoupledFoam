#!/usr/bin/env python3
"""Check T linear tolerances separately from normal solver termination."""
import argparse
import json
import math
import re
from pathlib import Path

PATTERN = re.compile(
    r'Solving for T, Initial residual = ([^,]+), Final residual = ([^,]+), No Iterations (\d+)'
)


def assess(text, tolerance=1e-8):
    solves = [dict(initial=float(a), final=float(b), iterations=int(n))
              for a, b, n in PATTERN.findall(text)]
    # Log rounding close to tolerance is not proof of convergence.
    failures = [i for i, s in enumerate(solves)
                if not math.isfinite(s['final']) or s['final'] > tolerance]
    phases = {}
    for line in text.splitlines():
        if line.startswith('M1_PHASE '):
            raw = dict(token.split('=', 1) for token in line.split()[1:])
            sample = {key: float(raw[key]) for key in ('time', 'massKg', 'alphaVolumeM3')}
            if 'nativeSensibleEnergyJ' in raw:
                sample['nativeSensibleEnergyJ'] = float(raw['nativeSensibleEnergyJ'])
            phases.setdefault(raw['phase'], []).append(sample)
    phase_inventory = {}
    for name, samples in phases.items():
        phase_inventory[name] = dict(samples=len(samples), initial=samples[0], final=samples[-1],
                                     mass_change_kg=samples[-1]['massKg']-samples[0]['massKg'],
                                     alpha_volume_change_m3=samples[-1]['alphaVolumeM3']-samples[0]['alphaVolumeM3'])
        if 'nativeSensibleEnergyJ' in samples[0] and 'nativeSensibleEnergyJ' in samples[-1]:
            phase_inventory[name]['native_sensible_energy_change_j'] = (
                samples[-1]['nativeSensibleEnergyJ']-samples[0]['nativeSensibleEnergyJ'])
    return dict(solves=solves, count=len(solves), failed_indices=failures,
                phase_inventory=phase_inventory,
                temperature_linear_check_passed=bool(solves) and not failures,
                normal_end=bool(re.search(r'^End\s*$', text, re.M)),
                tolerance=tolerance)


def closed_mass_check(text, relative_tolerance=1e-9):
    samples = []
    for line in text.splitlines():
        if line.startswith('M1_CONTINUUM '):
            raw = dict(token.split('=', 1) for token in line.split()[1:])
            samples.append({key: float(raw[key]) for key in
                            ('time', 'massKg', 'outwardMassFluxKgPerS')})
    if len(samples) < 2:
        raise ValueError('Closed mass check requires initial and final samples')
    flux_integral = 0.0
    maximum_residual = 0.0
    for before, after in zip(samples, samples[1:]):
        dt = after['time']-before['time']
        if dt <= 0:
            raise ValueError('Non-increasing diagnostic times')
        flux_integral += after['outwardMassFluxKgPerS']*dt
        maximum_residual = max(maximum_residual,
                               abs(after['massKg']-samples[0]['massKg']+flux_integral))
    delta = samples[-1]['massKg']-samples[0]['massKg']
    residual = delta+flux_integral
    relative = abs(residual)/abs(samples[0]['massKg'])
    maximum_relative = maximum_residual/abs(samples[0]['massKg'])
    return dict(initial_mass_kg=samples[0]['massKg'], final_mass_kg=samples[-1]['massKg'],
                inventory_change_kg=delta, right_endpoint_boundary_integral_kg=flux_integral,
                residual_kg=residual, relative_residual=relative,
                samples=len(samples), max_history_relative_residual=maximum_relative,
                relative_tolerance=relative_tolerance,
                closed_mass_check_passed=math.isfinite(maximum_relative) and maximum_relative <= relative_tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('--allow-upstream-reference-failure', action='store_true')
    args = parser.parse_args()
    cases = {}
    accepted = True
    for case in sorted(args.report.glob('T-*-average-*-subcycles-*')):
        log = case / 'log.compressibleLaserbeamFoam'
        text = log.read_text(errors='replace')
        result = assess(text)
        if (case/'phase-case.json').exists():
            result['closed_mass_balance'] = closed_mass_check(text)
            if not result['closed_mass_balance']['closed_mass_check_passed']:
                accepted = False
        reference = (args.allow_upstream_reference_failure
                     and case.name.startswith('T-upstream-'))
        result['reference_only'] = reference
        cases[case.name] = result
        if not result['normal_end'] or (not reference and not result['temperature_linear_check_passed']):
            accepted = False
        print(f"{case.name}: T solves={result['count']}, above tolerance={len(result['failed_indices'])}")
    accepted = accepted and bool(cases) and any(not c['reference_only'] for c in cases.values())
    summary = dict(schema_version=1, cases=cases, requested_gate_passed=accepted,
                   physical_validation='not_evaluated',
                   scope='T linear residuals and normal termination; closed synthetic phase cases also check mass. '
                         'Boundary integration is a right-endpoint estimate; energy balance not assessed.')
    (args.report / 'temperature-summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    if not accepted:
        raise SystemExit('Convergence or closed mass gate failed; report retained.')


if __name__ == '__main__':
    main()
