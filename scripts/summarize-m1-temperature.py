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
            phases.setdefault(raw['phase'], []).append(sample)
    phase_inventory = {}
    for name, samples in phases.items():
        phase_inventory[name] = dict(samples=len(samples), initial=samples[0], final=samples[-1],
                                     mass_change_kg=samples[-1]['massKg']-samples[0]['massKg'],
                                     alpha_volume_change_m3=samples[-1]['alphaVolumeM3']-samples[0]['alphaVolumeM3'])
    return dict(solves=solves, count=len(solves), failed_indices=failures,
                phase_inventory=phase_inventory,
                temperature_linear_check_passed=bool(solves) and not failures,
                normal_end=bool(re.search(r'^End\s*$', text, re.M)),
                tolerance=tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('--allow-upstream-reference-failure', action='store_true')
    args = parser.parse_args()
    cases = {}
    accepted = True
    for case in sorted(args.report.glob('T-*-average-*-subcycles-*')):
        log = case / 'log.compressibleLaserbeamFoam'
        result = assess(log.read_text(errors='replace'))
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
                   scope='Only T linear residuals and normal termination; other fields and physical energy balance not assessed.')
    (args.report / 'temperature-summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    if not accepted:
        raise SystemExit('Temperature convergence gate failed; report retained.')


if __name__ == '__main__':
    main()
