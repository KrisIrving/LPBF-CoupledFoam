#!/usr/bin/env python3
"""Summarize spatial integration gates; never apply the uniform energy gate."""
import argparse
import importlib.util
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('temperature', ROOT/'scripts/summarize-m1-temperature.py')
temperature = importlib.util.module_from_spec(spec)
spec.loader.exec_module(temperature)
spec_energy = importlib.util.spec_from_file_location('energy', ROOT/'scripts/summarize-m1-energy.py')
energy = importlib.util.module_from_spec(spec_energy)
spec_energy.loader.exec_module(energy)
FAMILIES = ('advection', 'spatial-evaporation', 'laser-heat', 'laser-phase')


def rows(text, prefix):
    return [dict(token.split('=',1) for token in line.split()[1:])
            for line in text.splitlines() if line.startswith(prefix+' ')]


def assess(case):
    meta = json.loads((case/'m1b-case.json').read_text())
    text = (case/'log.compressibleLaserbeamFoam').read_text(errors='replace')
    limits = meta['thresholds']
    linear = temperature.assess(text)
    continuum = rows(text,'M1_CONTINUUM')
    spatial = rows(text,'M1_SPATIAL')
    interfaces = rows(text,'M1_INTERFACE')
    expected_samples = round(meta['end_time_s']/meta['delta_t_s'])+1
    if len(continuum)!=expected_samples or len(interfaces)!=len(continuum) or len(spatial)!=9*len(continuum):
        raise ValueError('Missing or incomplete spatial observations')
    times = [float(r['time']) for r in continuum]
    if not math.isclose(times[0],0,abs_tol=1e-15) or any(
        not math.isclose(after-before,meta['delta_t_s'],rel_tol=1e-8,abs_tol=1e-15)
        for before,after in zip(times,times[1:])):
        raise ValueError('Missing or inconsistent physical time steps')
    if [float(r['time']) for r in interfaces]!=times:
        raise ValueError('Spatial and continuum times differ')
    phase_spatial = {}
    for row in spatial: phase_spatial.setdefault(row['phase'],[]).append(row)
    if any([float(r['time']) for r in series]!=times for series in phase_spatial.values()):
        raise ValueError('Phase spatial times differ')
    numeric = [(r,k) for r in continuum+spatial+interfaces for k in r if k not in ('phase','solver')]
    if any(not math.isfinite(float(r[k])) for r,k in numeric):
        raise ValueError('Non-finite diagnostics')
    mass = temperature.closed_mass_check(text,limits['closed_mass_relative_tolerance'])
    minimum = min(float(r['alphaMin']) for r in spatial)
    maximum = max(float(r['alphaMax']) for r in spatial)
    sum_error = max(max(abs(float(r['alphaSumMin'])-1),abs(float(r['alphaSumMax'])-1)) for r in interfaces)
    normal = linear['normal_end'] and not re.search(r'FOAM FATAL|SIGFPE signal|(^|[^a-z])(nan|inf)([^a-z]|$)',text,re.I)
    end_reached = math.isclose(times[-1],meta['end_time_s'],rel_tol=1e-9,abs_tol=1e-15)
    gates = dict(normal_end=bool(normal),end_time_reached=end_reached,
                 temperature_linear=linear['temperature_linear_check_passed'],
                 closed_mass=mass['closed_mass_check_passed'],
                 alpha_bounds=minimum>=-limits['alpha_bound_tolerance'] and maximum<=1+limits['alpha_bound_tolerance'],
                 alpha_sum=sum_error<=limits['alpha_sum_tolerance'])
    absorbed = math.fsum(float(after['absorbedLaserPowerW'])*(float(after['time'])-float(before['time']))
                         for before,after in zip(continuum,continuum[1:]))
    extra = {}
    optics=rows(text,'M1_OPTICS')
    if optics:
        if [float(r['time']) for r in optics]!=times:
            raise ValueError('Optical diagnostics and continuum times differ')
        resistivities=[float(r['condensedWeightedResistivityOhmM']) for r in optics]
        if any(not math.isfinite(x) for x in resistivities):
            raise ValueError('Non-finite optical diagnostic')
        extra['final_condensed_weighted_resistivity_ohm_m']=resistivities[-1]
    elif meta.get('condensed_phase_optics',False):
        raise ValueError('Missing optical diagnostics for condensed-phase candidate')
    if meta['family']=='advection':
        series = phase_spatial['metal1']
        error = float(interfaces[-1]['analyticMeanAbsError'])
        initial_volume = float(series[0]['alphaVolumeM3'])
        volume_change = abs(float(series[-1]['alphaVolumeM3'])-initial_volume)/initial_volume
        gates.update(analytic_advection=0<=error<=limits['advection_mean_abs_error'],
                     phase_volume=volume_change<=limits['advection_phase_volume_relative_tolerance'])
        extra.update(final_analytic_mean_abs_error=error,phase_volume_relative_change=volume_change)
    if meta['laser_power_w']>0:
        gates['positive_absorbed_energy']=absorbed>0
        gates['temperature_response']=max(float(r['TmaxK']) for r in continuum)>meta['initial_temperature_K']+1e-8
    phase_rows = {}
    for row in rows(text,'M1_PHASE'): phase_rows.setdefault(row['phase'],[]).append(row)
    if set(phase_rows)!=set(phase_spatial) or any([float(r['time']) for r in series]!=times for series in phase_rows.values()):
        raise ValueError('Phase inventory times differ')
    offset=meta['vapour_reference_energy_offset_J_per_kg']
    energies=[math.fsum(float(series[i]['nativeSensibleEnergyJ']) for series in phase_rows.values())
              +offset*float(phase_rows['metal1vapour'][i]['massKg'])+float(continuum[i]['kineticEnergyJ'])
              for i in range(len(times))]
    if any(not math.isfinite(e) for e in energies):
        raise ValueError('Non-finite phase energy observations')
    mv=phase_rows['metal1vapour']
    delta_vapour=float(mv[-1]['massKg'])-float(mv[0]['massKg'])
    phase_change_source = [float(r['temperatureSourceKgKPerS']) for r in rows(text,'M1_COUPLING') if r['stage']=='afterTemperature']
    if not phase_change_source or any(not math.isfinite(x) for x in phase_change_source):
        raise ValueError('Missing/non-finite phase source observations')
    if not meta['phase_change_enabled']:
        gates['disabled_phase_source_zero']=max(abs(x) for x in phase_change_source)<=1e-20
    else:
        gates['active_phase_source_observed']=max(abs(x) for x in phase_change_source)>limits['minimum_active_temperature_source_kg_k_per_s']
    return dict(metadata=meta,gates=gates,passed=all(gates.values()),
                temperature_equation_budget=energy.temperature_budget(text,continuum,meta,offset,delta_vapour)
                    if meta.get('common_constant_cv_J_per_kg_K') else dict(available=False,reason='Missing explicit Cv metadata'),
                closed_mass_balance=mass,alpha_min=minimum,alpha_max=maximum,alpha_sum_max_error=sum_error,
                observations=dict(absorbed_laser_energy_j=absorbed,
                    configured_incident_energy_j=meta['laser_power_w']*(times[-1]-times[0]),
                    reference_energy_change_j=energies[-1]-energies[0],
                    energy_change_minus_absorbed_j=energies[-1]-energies[0]-absorbed,
                    vapour_mass_change_kg=delta_vapour,
                    max_abs_temperature_source_kg_k_per_s=max(abs(x) for x in phase_change_source),
                    final_Tmin_K=float(continuum[-1]['TminK']),final_Tmax_K=float(continuum[-1]['TmaxK']),
                    **extra),physical_validation='not_evaluated')


def summarize(report):
    cases={}
    for family in FAMILIES:
        for level in ('coarse','fine'):
            name=family+'-'+level
            try:
                status=(report/name/'case-result.txt').read_text().strip()
                if status!='exit_status=0':
                    detail=''
                    log=report/name/'log.compressibleLaserbeamFoam'
                    if log.exists():
                        match=re.search(r"Entry '[^']+' not found[^\n]*",log.read_text(errors='replace'))
                        if match:detail='; '+match.group()
                    raise ValueError('Case execution failed: '+status+detail)
                cases[name]=assess(report/name)
            except (OSError,ValueError,KeyError,ZeroDivisionError,IndexError) as error:
                cases[name]=dict(passed=False,error=str(error))
    comparisons={}
    for family in FAMILIES:
        coarse,fine=cases[family+'-coarse'],cases[family+'-fine']
        if 'observations' in coarse and 'observations' in fine:
            comparisons[family]={key:dict(coarse=coarse['observations'][key],fine=fine['observations'][key],
                                        fine_minus_coarse=fine['observations'][key]-coarse['observations'][key])
                                 for key in coarse['observations']}
    return dict(schema_version=1,cases=cases,resolution_comparisons=comparisons,
                requested_gate_passed=all(c['passed'] for c in cases.values()),
                physical_validation='not_evaluated',
                scope='Spatial synthetic integration only. Energy minus absorbed laser energy is an '
                      'observation, not a conservation gate: pressure, phase/fusion transport and boundary '
                      'energy accounting are not yet complete. Two resolutions do not prove convergence. '
                      'Nonuniform phase-source tests exercise upstream volume-distributed kinetics, '
                      'not a validated interface-localized evaporation law.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    args=parser.parse_args()
    result=summarize(args.report)
    (args.report/'m1b-summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for name,c in result['cases'].items(): print(name, 'PASS' if c['passed'] else 'FAIL',c.get('error',''))
    if not result['requested_gate_passed']:raise SystemExit('Integrated gates failed; all case feedback retained.')
