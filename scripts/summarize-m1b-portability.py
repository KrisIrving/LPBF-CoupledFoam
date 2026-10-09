#!/usr/bin/env python3
"""Compare whole diagnostic histories; thresholds fixed before user execution."""
import argparse
import importlib.util
import json
import math
import re
from pathlib import Path

spec = importlib.util.spec_from_file_location('integrated', Path(__file__).with_name('summarize-m1b.py'))
integrated = importlib.util.module_from_spec(spec)
spec.loader.exec_module(integrated)
LIMITS = dict(mass_relative=1e-8, phase_volume_domain_relative=1e-8,
              temperature_absolute_K=1e-4, absorbed_absolute_J=1e-14,
              energy_change_absolute_J=1e-13, energy_signal_relative=1e-5)


def merge(first, second):
    a, b = first.read_text(errors='replace'), second.read_text(errors='replace')
    for text in (a, b):
        if not re.search(r'^End$', text, re.M):
            raise ValueError('Both restart segments must end normally')
    ar, br = integrated.rows(a,'M1_CONTINUUM'), integrated.rows(b,'M1_CONTINUUM')
    midpoint = float(ar[-1]['time'])
    if float(br[0]['time']) != midpoint or float(br[1]['time']) <= midpoint:
        raise ValueError('Restart times do not join at the checkpoint')
    prefixes = ('M1_CONTINUUM','M1_PHASE','M1_SPATIAL','M1_INTERFACE','M1_OPTICS')
    kept = []
    for line in b.splitlines():
        if any(line.startswith(p+' ') for p in prefixes):
            tokens = dict(t.split('=',1) for t in line.split()[1:])
            if float(tokens['time']) == midpoint:
                continue
        kept.append(line)
    return re.sub(r'^End$', 'M1_SEGMENT_END', a, flags=re.M)+'\n'+'\n'.join(kept)+'\n'


def history(case):
    text = (case/'log.compressibleLaserbeamFoam').read_text(errors='replace')
    meta = json.loads((case/'m1b-case.json').read_text())
    continuum = integrated.rows(text,'M1_CONTINUUM')
    phases = {}
    for r in integrated.rows(text,'M1_PHASE'): phases.setdefault(r['phase'],[]).append(r)
    spatial = {}
    for r in integrated.rows(text,'M1_SPATIAL'): spatial.setdefault(r['phase'],[]).append(r)
    cumulative = 0.0
    result = []
    for i, r in enumerate(continuum):
        if i:
            cumulative += float(r['absorbedLaserPowerW'])*(float(r['time'])-float(continuum[i-1]['time']))
        native = math.fsum(float(s[i]['nativeSensibleEnergyJ']) for s in phases.values())
        native += meta['vapour_reference_energy_offset_J_per_kg']*float(phases['metal1vapour'][i]['massKg'])
        native += float(r['kineticEnergyJ'])
        result.append(dict(time=float(r['time']),mass=float(r['massKg']),
                           Tmin=float(r['TminK']),Tmax=float(r['TmaxK']),absorbed=cumulative,energy=native,
                           volumes={p:float(s[i]['alphaVolumeM3']) for p,s in spatial.items()}))
    return meta,result


def compare(base, candidate):
    meta,a = history(base)
    other,b = history(candidate)
    if meta != other or [r['time'] for r in a] != [r['time'] for r in b]:
        raise ValueError('Physical configuration or sample times differ')
    differences = dict(mass_relative=max(abs(x['mass']-y['mass'])/a[0]['mass'] for x,y in zip(a,b)),
        phase_volume_domain_relative=max(abs(x['volumes'][p]-y['volumes'][p])/meta['length_m']**3
                                        for x,y in zip(a,b) for p in x['volumes']),
        temperature_absolute_K=max(abs(x[k]-y[k]) for x,y in zip(a,b) for k in ('Tmin','Tmax')),
        absorbed_absolute_J=max(abs(x['absorbed']-y['absorbed']) for x,y in zip(a,b)),
        energy_change_absolute_J=max(abs((x['energy']-a[0]['energy'])-(y['energy']-b[0]['energy']))
                                     for x,y in zip(a,b)))
    gates = {k:v <= LIMITS[k]+(LIMITS['energy_signal_relative']*abs(a[-1]['absorbed'])
                              if k in ('absorbed_absolute_J','energy_change_absolute_J') else 0)
             for k,v in differences.items()}
    return dict(max_history_differences=differences,gates=gates,passed=all(gates.values()))


def summarize(report):
    cases, comparisons = {}, {}
    for name in ('continuous','restart','mpi2'):
        try:
            if (report/name/'case-result.txt').read_text().strip() != 'exit_status=0':
                raise ValueError('Case execution failed')
            cases[name] = integrated.assess(report/name)
        except (OSError,ValueError,KeyError,IndexError,ZeroDivisionError) as error:
            cases[name] = dict(passed=False,error=str(error))
    for name in ('restart','mpi2'):
        if cases['continuous']['passed'] and cases[name]['passed']:
            try: comparisons[name] = compare(report/'continuous',report/name)
            except (OSError,ValueError,KeyError,IndexError,ZeroDivisionError) as error:
                comparisons[name] = dict(passed=False,error=str(error))
        else: comparisons[name] = dict(passed=False,error='Prerequisite case gates failed; comparison not evaluated')
    return dict(package='M1B-02',thresholds=LIMITS,cases=cases,comparisons=comparisons,
                requested_gate_passed=all(c['passed'] for c in list(cases.values())+list(comparisons.values())),
                physical_validation='not_evaluated',
                scope='Synthetic weak-phase laser fixture only; global inventories and temperature extrema. '
                      'No fieldwise equivalence, strong evaporation, full energy closure or material validation.')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    parser.add_argument('--merge',action='store_true')
    args = parser.parse_args()
    if args.merge:
        (args.report/'log.compressibleLaserbeamFoam').write_text(merge(args.report/'log.first',args.report/'log.restart'))
    else:
        result = summarize(args.report)
        (args.report/'m1b-portability-summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        for name,c in {**result['cases'],**{'compare-'+k:v for k,v in result['comparisons'].items()}}.items():
            print(name,'PASS' if c['passed'] else 'FAIL',c.get('error',''))
        if not result['requested_gate_passed']: raise SystemExit('Portability gates failed; feedback retained.')
