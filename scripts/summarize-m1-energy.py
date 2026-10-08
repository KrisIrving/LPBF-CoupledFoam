#!/usr/bin/env python3
"""Observe closed energy with a fixed phase offset calibrated at equilibrium."""
import argparse
import json
import math
from pathlib import Path


def rows(text, prefix):
    result = []
    for line in text.splitlines():
        if line.startswith(prefix+' '):
            raw = dict(token.split('=', 1) for token in line.split()[1:])
            result.append(raw)
    return result


def temperature_budget(text, continuum, meta, offset, delta_mv):
    cv = meta.get('common_constant_cv_J_per_kg_K')
    budget_rows = rows(text, 'M1_T_BUDGET')
    if not budget_rows:
        return dict(available=False)
    if meta.get('thermo_fixture') != 'consistent' or cv is None or not math.isfinite(cv) or cv <= 0:
        raise ValueError('Temperature-budget energy conversion requires common constant Cv metadata')
    # PIMPLE repeatedly solves the same time interval. Only its final thermal
    # correction contributes once to this diagnostic integral.
    final_by_time = {float(r['time']): r for r in budget_rows}
    times = [float(r['time']) for r in continuum]
    if set(final_by_time) != set(times[1:]):
        raise ValueError('Temperature budget and continuum times are not aligned')
    keys = ('storage', 'transport', 'diffusion', 'mechanical', 'laser', 'fusion', 'phaseChange', 'residual')
    integrals = {k: 0.0 for k in keys}
    history = []
    for before, after in zip(times, times[1:]):
        if after <= before:
            raise ValueError('Non-increasing temperature budget times')
        sample = final_by_time[after]
        for key in keys:
            value = float(sample[key+'KgKPerS'])
            if not math.isfinite(value):
                raise ValueError('Non-finite temperature budget term')
            integrals[key] += cv*value*(after-before)
        history.append(dict(time_s=after, inventory_kg_k=float(sample['inventoryKgK']),
                            old_inventory_kg_k=float(sample['oldInventoryKgK']),
                            terms_kg_k_per_s={k: float(sample[k+'KgKPerS']) for k in keys}))
    return dict(available=True, constant_cv_J_per_kg_K=cv,
                final_outer_temperature_samples=history,
                integrated_signed_equation_terms_j=integrals,
                reference_offset_inventory_change_j=offset*delta_mv,
                phase_source_inventory_mismatch_j=offset*delta_mv-integrals['phaseChange'],
                scope='Unrelaxed rho*T equation terms evaluated before mixture.correct and pressure '
                      'correction. Last thermal solve per time is integrated once. Constant Cv conversion '
                      'is restricted to this fixture; EOS native energy after pressure is a separate ledger. '
                      'Reconstructed explicit transport can differ from implicit matrix transport in '
                      'nonuniform cases; this is a diagnostic, not a new conservation gate.')


def audit(report, gate=False, absolute_tolerance=1e-9, relative_tolerance=1e-4):
    if any(not math.isfinite(x) or x < 0 for x in (absolute_tolerance, relative_tolerance)):
        raise ValueError('Energy tolerances must be finite and nonnegative')
    cases = []
    for case in sorted(report.glob('T-*-average-*-subcycles-*')):
        meta = case/'phase-case.json'
        if meta.exists():
            log = (case/'log.compressibleLaserbeamFoam').read_text()
            phases = {}
            for row in rows(log, 'M1_PHASE'):
                phases.setdefault(row['phase'], []).append(row)
            cases.append((case.name, json.loads(meta.read_text()), phases,
                          rows(log, 'M1_CONTINUUM')))
    refs = [c for c in cases if c[1]['direction'] == 'equilibrium']
    if not refs:
        raise ValueError('Equilibrium reference required')
    _, reference, phases, _ = refs[0]
    l, v = phases['metal1'][0], phases['metal1vapour'][0]
    el = float(l['nativeSensibleEnergyJ'])/float(l['massKg'])
    ev = float(v['nativeSensibleEnergyJ'])/float(v['massKg'])
    specific_volume_l = float(l['alphaVolumeM3'])/float(l['massKg'])
    specific_volume_v = float(v['alphaVolumeM3'])/float(v['massKg'])
    p = reference['initial_pressure_Pa']
    latent = reference['latent_heat_J_per_kg']
    # Lv is an enthalpy gap: h_v-h_l = e_v-e_l + p*(v_v-v_l).
    offset = latent-(ev-el)-p*(specific_volume_v-specific_volume_l)
    for _, ref_meta, ref_phases, _ in refs[1:]:
        for phase in ('metal1', 'metal1vapour'):
            for key in ('massKg', 'alphaVolumeM3', 'nativeSensibleEnergyJ'):
                if not math.isclose(float(ref_phases[phase][0][key]),
                                    float(phases[phase][0][key]), rel_tol=1e-12):
                    raise ValueError('Equilibrium initial references differ')
        if ref_meta['latent_heat_J_per_kg'] != latent or ref_meta['initial_pressure_Pa'] != p:
            raise ValueError('Equilibrium latent/pressure references differ')
    assert math.isclose(ev+offset-el+p*(specific_volume_v-specific_volume_l),
                        latent, rel_tol=1e-8, abs_tol=1e-8)
    results = {}
    for name, meta, inventories, continuum in cases:
        configured_offset = meta.get('vapour_reference_energy_offset_J_per_kg')
        if meta.get('referenceInternalEnergySource', False) and (
            configured_offset is None or not math.isclose(configured_offset, offset, rel_tol=1e-10)
        ):
            raise ValueError('Configured source offset differs from equilibrium ledger calibration')
        energies = []
        for index, sample in enumerate(continuum):
            time = float(sample['time'])
            if any(len(series) != len(continuum) or float(series[index]['time']) != time
                   for series in inventories.values()):
                raise ValueError('Phase and continuum samples are not aligned')
            native = math.fsum(float(series[index]['nativeSensibleEnergyJ'])
                               for series in inventories.values())
            vapour_mass = float(inventories['metal1vapour'][index]['massKg'])
            energy = native+offset*vapour_mass+float(sample['kineticEnergyJ'])
            energies.append((time, energy))
        mv = inventories['metal1vapour']
        delta_mv = float(mv[-1]['massKg'])-float(mv[0]['massKg'])
        history = []
        for index, (time, energy) in enumerate(energies):
            change = energy-energies[0][1]
            transfer = abs(latent*(float(mv[index]['massKg'])-float(mv[0]['massKg'])))
            tolerance = absolute_tolerance+relative_tolerance*transfer
            history.append(dict(time_s=time, reference_energy_change_j=change,
                                latent_transfer_scale_j=transfer, tolerance_j=tolerance,
                                within_tolerance=math.isfinite(change) and abs(change) <= tolerance))
        results[name] = dict(direction=meta['direction'],
            thermo_fixture=meta.get('thermo_fixture', 'legacy'),
            reference_internal_energy_source=meta.get('referenceInternalEnergySource', False),
            common_latent_heat_source=meta.get('commonLatentHeatSource', False), samples=len(energies),
            temperature_equation_budget=temperature_budget(
                (report/name/'log.compressibleLaserbeamFoam').read_text(), continuum, meta, offset, delta_mv),
            initial_reference_energy_j=energies[0][1], final_reference_energy_j=energies[-1][1],
            reference_energy_change_j=energies[-1][1]-energies[0][1],
            max_abs_reference_energy_change_j=max(abs(e-energies[0][1]) for _, e in energies),
            vapour_mass_change_kg=delta_mv, latent_transfer_scale_j=abs(latent*delta_mv),
            energy_history=history,
            energy_regression_passed=all(s['within_tolerance'] for s in history) if gate else None)
    return dict(schema_version=2, vapour_constant_energy_offset_J_per_kg=offset,
                reference_temperature_K=reference['temperature_K'], reference_pressure_Pa=p,
                cases=results, physical_energy_validation='not_evaluated',
                energy_regression_gate=dict(enabled=gate, absolute_tolerance_j=absolute_tolerance,
                    relative_transfer_tolerance=relative_tolerance,
                    passed=all(c['energy_regression_passed'] for c in results.values()) if gate else None),
                scope='Synthetic closed adiabatic zero-laser case only. Offset calibrated once from '
                      'equilibrium native thermo energies so the reference enthalpy gap equals Lv. '
                      'The ledger does not modify solver states; case source modes are recorded separately. '
                      'Optional regression gate checks every logged sample against an absolute floor '
                      'plus a tolerance relative to net Lv*delta(vapour mass), not initial total energy. '
                      'This scale is diagnostic, not the internal-energy source. '
                      'It is intended for monotonic uniform synthetic phase exchange only. '
                      'EOS thermodynamic consistency and off-reference latent behavior remain to be audited.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    parser.add_argument('--gate', action='store_true', help='Require the synthetic energy regression gate')
    parser.add_argument('--absolute-tolerance', type=float, default=1e-9, help='Printed-energy noise floor in J')
    parser.add_argument('--relative-tolerance', type=float, default=1e-4, help='Fraction of net latent transfer scale')
    args = parser.parse_args()
    result = audit(args.report, args.gate, args.absolute_tolerance, args.relative_tolerance)
    (args.report/'energy-summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    if args.gate and not result['energy_regression_gate']['passed']:
        raise SystemExit('Synthetic energy regression failed; see energy-summary.json.')
    print('Reference-calibrated energy summary written; physical validation pending.')
