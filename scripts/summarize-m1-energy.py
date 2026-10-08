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


def audit(report):
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
    if len(refs) != 1:
        raise ValueError('Exactly one equilibrium reference required')
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
    assert math.isclose(ev+offset-el+p*(specific_volume_v-specific_volume_l),
                        latent, rel_tol=1e-8, abs_tol=1e-8)
    results = {}
    for name, meta, inventories, continuum in cases:
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
        results[name] = dict(direction=meta['direction'], samples=len(energies),
            initial_reference_energy_j=energies[0][1], final_reference_energy_j=energies[-1][1],
            reference_energy_change_j=energies[-1][1]-energies[0][1],
            max_abs_reference_energy_change_j=max(abs(e-energies[0][1]) for _, e in energies),
            vapour_mass_change_kg=delta_mv, latent_transfer_scale_j=abs(latent*delta_mv))
    return dict(schema_version=1, vapour_constant_energy_offset_J_per_kg=offset,
                reference_temperature_K=reference['temperature_K'], reference_pressure_Pa=p,
                cases=results, physical_energy_validation='not_evaluated',
                scope='Synthetic closed adiabatic zero-laser case only. Offset calibrated once from '
                      'equilibrium native thermo energies so the reference enthalpy gap equals Lv. '
                      'Solver energies and sources are unchanged; no energy tolerance gate is imposed. '
                      'EOS thermodynamic consistency and off-reference latent behavior remain to be audited.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    args = parser.parse_args()
    result = audit(args.report)
    (args.report/'energy-summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Reference-calibrated energy observations written; physical validation pending.')
