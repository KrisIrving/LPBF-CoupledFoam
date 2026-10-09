#!/usr/bin/env python3
"""Validate paired interface sources, inventory/boundary balances and response."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

spec=importlib.util.spec_from_file_location('integrated', Path(__file__).with_name('summarize-m1b.py'))
integrated=importlib.util.module_from_spec(spec)
spec.loader.exec_module(integrated)
CASES=('off','prescribed-coarse','prescribed-fine','prescribed-half-dt','kinetic','condensation','open')


def last_rows(text, prefix, key='time'):
    return {row[key]:row for row in integrated.rows(text,prefix)}


def assess(case):
    if (case/'case-result.txt').read_text().strip()!='exit_status=0':
        raise ValueError('Case execution failed')
    meta=json.loads((case/'m1c-case.json').read_text())
    text=(case/'log.compressibleLaserbeamFoam').read_text(errors='replace')
    limits=meta['thresholds']
    continuum=integrated.rows(text,'M1_CONTINUUM')
    if len(continuum)!=round(meta['end_time_s']/meta['delta_t_s'])+1:
        raise ValueError('Incomplete continuum history')
    times=[float(row['time']) for row in continuum]
    if times[0]!=0 or not math.isclose(times[-1],meta['end_time_s'],rel_tol=1e-9) or any(
        not math.isclose(b-a,meta['delta_t_s'],rel_tol=1e-8,abs_tol=1e-18)
        for a,b in zip(times,times[1:])):
        raise ValueError('Missing/inconsistent time steps')
    sources=last_rows(text,'M1C_SOURCE')
    flows=last_rows(text,'M1C_FLOW')
    interfaces=last_rows(text,'M1_INTERFACE')
    phases={}
    for row in integrated.rows(text,'M1_PHASE'):
        phases.setdefault(row['phase'],{})[row['time']]=row
    fluxes={}
    for row in integrated.rows(text,'M1C_PHASE_FLUX'):
        fluxes.setdefault(row['phase'],{})[row['time']]=row
    keys=[row['time'] for row in continuum]
    steps=set(keys[1:])
    spatial=integrated.rows(text,'M1_SPATIAL')
    spatial_phases={}
    for row in spatial: spatial_phases.setdefault(row['phase'],{})[row['time']]=row
    if (set(sources)!=steps or set(flows)!=set(keys) or set(interfaces)!=set(keys)
        or len(phases)!=9 or set(phases)!=set(fluxes) or set(phases)!=set(spatial_phases)
        or any(set(rows)!=set(keys) for rows in phases.values())
        or any(set(rows)!=set(keys) for rows in spatial_phases.values())
        or any(set(rows)!=steps for rows in fluxes.values())):
        raise ValueError('Incomplete source, phase, boundary or pressure diagnostics')
    numeric=continuum+list(sources.values())+list(flows.values())+list(interfaces.values())+spatial
    numeric += [row for series in phases.values() for row in series.values()]
    numeric += [row for series in fluxes.values() for row in series.values()]
    if any(not math.isfinite(float(v)) for row in numeric for k,v in row.items()
           if k not in ('solver','phase','model','enabled')):
        raise ValueError('Non-finite diagnostic')
    if any(row['model']!=meta['model'] or (row['enabled'].lower() in ('1','true','on'))!=meta['enabled']
           or not math.isclose(float(row['dt']),meta['delta_t_s'],rel_tol=1e-10)
           for row in sources.values()):
        raise ValueError('Source control/metadata mismatch')
    source_rows=[sources[k] for k in keys[1:]]
    # Each outer iteration replaces a trial step; integrate its final source once.
    transferred=boundary_total=boundary_liquid=boundary_vapour=0.0
    residuals={'total':[], 'liquid':[], 'vapour':[]}
    initial_mass=float(continuum[0]['massKg'])
    initial_liquid=float(phases['metal1'][keys[0]]['massKg'])
    initial_vapour=float(phases['metal1vapour'][keys[0]]['massKg'])
    transfer_history=[]
    for index,key in enumerate(keys[1:],1):
        dt=times[index]-times[index-1]
        transferred+=float(sources[key]['transferredKgPerS'])*dt
        boundary_total+=float(continuum[index]['outwardMassFluxKgPerS'])*dt
        boundary_liquid+=float(fluxes['metal1'][key]['outwardKgPerS'])*dt
        boundary_vapour+=float(fluxes['metal1vapour'][key]['outwardKgPerS'])*dt
        residuals['total'].append(float(continuum[index]['massKg'])-initial_mass+boundary_total)
        residuals['liquid'].append(float(phases['metal1'][key]['massKg'])-initial_liquid+boundary_liquid+transferred)
        residuals['vapour'].append(float(phases['metal1vapour'][key]['massKg'])-initial_vapour+boundary_vapour-transferred)
        transfer_history.append(transferred)
    budget=limits['mass_balance_absolute_kg']+limits['mass_balance_signal_fraction']*meta['comparison_mass_scale_kg']
    max_residuals={name:max(map(abs,values)) for name,values in residuals.items()}
    area=meta['expected_planar_area_m2']
    first_source=integrated.rows(text,'M1C_SOURCE')[0]  # Before any first-step phase/thermal update.
    initial_area_error=abs(float(first_source['areaM2'])/area-1)
    area_error=max(abs(float(row['areaM2'])/area-1) for row in source_rows)
    pair_error=max(abs(float(row['pairMassResidualKgPerS'])) for row in source_rows)
    requested=max(abs(float(row['requestedKgPerS'])) for row in source_rows)
    limiting=max(float(row['limitedAbsKgPerS']) for row in source_rows)
    energy_error=max(abs(float(row['phasePowerW'])-float(row['transferredKgPerS'])*meta['energy_jump_J_kg'])
                     for row in source_rows)
    rate_scale=meta['comparison_mass_scale_kg']/meta['end_time_s']
    expected_j=meta['prescribed_mass_flux_kg_m2_s']
    if meta['model']=='kinetic':
        expected_j=meta['accommodation']*(meta['reference_pressure_Pa']-meta['initial_pressure_Pa'])*math.sqrt(
            meta['molar_mass_kg_mol']/(2*math.pi*meta['gas_constant_J_mol_K']*meta['initial_temperature_K']))
    expected_initial=expected_j*float(first_source['areaM2']) if meta['enabled'] else 0
    initial_closure_error=abs(float(first_source['requestedKgPerS'])-expected_initial)
    prescribed_error=max(abs(float(row['requestedKgPerS'])-(
        meta['prescribed_mass_flux_kg_m2_s']*float(row['areaM2']) if meta['enabled'] else 0)) for row in source_rows)
    linear=integrated.temperature.assess(text)
    pressure_rise=float(flows[keys[-1]]['pressureMeanPa'])-float(flows[keys[0]]['pressureMeanPa'])
    speed=max(float(row['maxSpeedMS']) for row in continuum)
    alpha_min=min(float(row['alphaMin']) for row in spatial)
    alpha_max=max(float(row['alphaMax']) for row in spatial)
    alpha_sum=max(max(abs(float(row['alphaSumMin'])-1),abs(float(row['alphaSumMax'])-1)) for row in interfaces.values())
    gates=dict(normal_end=linear['normal_end'] and 'FOAM FATAL' not in text,
               temperature_linear=linear['temperature_linear_check_passed'],
               alpha_bounds=alpha_min>=-limits['alpha_bound'] and alpha_max<=1+limits['alpha_bound'],
               alpha_sum=alpha_sum<=limits['alpha_sum'],
               initial_area=initial_area_error<=limits['initial_area_relative'],
               area_history=area_error<=limits['max_area_relative'],
               initial_closure=initial_closure_error<=rate_scale*1e-10,
               paired_mass_source=pair_error<=rate_scale*limits['source_pair_relative'],
               shared_energy_source=energy_error<=rate_scale*abs(meta['energy_jump_J_kg'])*limits['source_energy_relative'],
               limiter_inactive=limiting<=max(requested,rate_scale)*limits['limiter_fraction'],
               **{name+'_inventory_balance':value<=budget for name,value in max_residuals.items()})
    if meta['model']=='prescribed': gates['prescribed_flux']=prescribed_error<=rate_scale*1e-10
    if meta['enabled']:
        direction=-1 if meta['initial_pressure_Pa']>meta['reference_pressure_Pa'] and meta['model']=='kinetic' else 1
        gates['transfer_direction']=direction*transferred>0
        gates['gas_velocity_response']=speed>=limits['minimum_speed_m_s']
        if not meta['open_boundary']: gates['pressure_response']=direction*pressure_rise>=limits['required_pressure_rise_Pa']
        else: gates['outlet_mass_response']=boundary_vapour>0
    else: gates['disabled_source']=transferred==0 and requested==0 and energy_error==0
    return dict(metadata=meta,gates=gates,passed=all(gates.values()),
                observations=dict(transferred_mass_kg=transferred,gas_outward_mass_kg=boundary_vapour,
                    initial_area_relative_error=initial_area_error,max_area_relative_error=area_error,
                    max_inventory_residual_kg=max_residuals,mass_balance_budget_kg=budget,
                    max_pair_source_residual_kg_s=pair_error,max_energy_source_error_w=energy_error,
                    pressure_mean_rise_pa=pressure_rise,max_speed_m_s=speed),
                physical_validation='not_evaluated')


def summarize(report):
    cases={}
    for name in CASES:
        try: cases[name]=assess(report/name)
        except (OSError,ValueError,KeyError,IndexError,ZeroDivisionError) as error:
            cases[name]=dict(passed=False,error=str(error))
    # Resolution differences are observations, never inferred convergence order.
    sensitivity={}
    baseline=cases['prescribed-coarse']
    if baseline['passed']:
        for name in ('prescribed-fine','prescribed-half-dt'):
            if cases[name]['passed']:
                sensitivity[name]={key:cases[name]['observations'][key]-baseline['observations'][key]
                    for key in ('transferred_mass_kg','pressure_mean_rise_pa','max_speed_m_s')}
    return dict(package='M1C-01',revision=1,cases=cases,resolution_differences=sensitivity,
                requested_gate_passed=all(case['passed'] for case in cases.values()),
                physical_validation='not_evaluated',
                scope='Paired interface mass/energy and gas response on synthetic pure vapour; '
                      'no explicit exit momentum, shielding gas, material or resolved-particle validation.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    args=parser.parse_args()
    result=summarize(args.report)
    (args.report/'m1c-interface-summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for name,case in result['cases'].items():
        failed=', '.join(k for k,v in case.get('gates',{}).items() if not v)
        print(name,'PASS' if case['passed'] else 'FAIL',case.get('error',failed))
    if not result['requested_gate_passed']: raise SystemExit('Interface gates failed; feedback retained.')
