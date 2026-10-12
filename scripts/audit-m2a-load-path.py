"""Offline separation of wall constraints, independent sampling and load ledger.

Reads retained trajectories only. Does not run CFD or replace physical gates.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def audit_case(case):
    meta=json.loads((case/'M2A_META.json').read_text())
    with (case/'mechanical-history.csv').open(newline='') as stream:
        rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]
    if not rows or any(not math.isfinite(v) for row in rows for v in row.values()):
        raise ValueError('Missing or nonfinite trajectory')
    complete=(case/'result.txt').read_text().strip()=='exit_status=0' and len(rows)==meta['expected_rows']
    mode=meta['mode'];radius=meta['radius'];fixed=mode=='fixed'
    if mode not in ('fixed','rotate'):raise ValueError('Static sphere audit only')
    prefix='f' if fixed else 't';axis='x' if fixed else 'z';sign=1 if fixed else -1
    reference=(6*math.pi*meta['rho']*meta['nu']*radius*meta['reference_flow_speed'] if fixed
               else 8*math.pi*meta['rho']*meta['nu']*radius**3*meta['reference_rotation_speed'])
    speed=meta['reference_flow_speed'] if fixed else radius*meta['reference_rotation_speed']
    tail=rows[-10:]
    mean=lambda key:sum(row[key] for row in tail)/len(tail)
    ratios={key:sign*mean(column)/reference for key,column in (
        ('constraint_reaction','constraint_'+prefix+axis),('inventory_correction','inertia_'+prefix+axis),
        ('ledger',prefix+axis),('pressure_traction','stress_pressure_'+prefix+axis),
        ('viscous_traction','stress_viscous_'+prefix+axis))}
    ratios['constraint_reaction']*=-1
    ratios['traction']=ratios['pressure_traction']+ratios['viscous_traction']
    identity=max(abs(row[prefix+c]+row['constraint_'+prefix+c]-row['inertia_'+prefix+c])
                 for row in rows for c in 'xyz')/reference
    return {'completed_windows':len(rows),'complete':complete,
        'sample_scope':'last10 complete windows' if complete else 'last up-to10 retained windows; diagnostic only',
        'tail_windows':len(tail),'tail_start_s':tail[0]['time'],'tail_end_s':tail[-1]['time'],
        'normalised_load_components':ratios,'max_load_identity_defect':identity,
        'max_collocation_residual_m_s':max(row['joint_wall_residual'] for row in rows),
        'max_independent_surface_slip_relative':max(row['slip_rms'] for row in rows)/speed,
        'max_bulk_slip_relative':max(row['bulk_target_slip_rms'] for row in rows)/speed,
        'physical_pass':'not_assessed_by_this_diagnostic'}


def run(report):
    result={'scope':__doc__,'cases':{}}
    for meta in sorted(report.glob('*/M2A_META.json')):
        try:entry=audit_case(meta.parent)
        except (OSError,ValueError,KeyError) as error:entry={'error':str(error),'physical_pass':'not_assessed_by_this_diagnostic'}
        result['cases'][meta.parent.name]=entry
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(run(args.report),indent=2,allow_nan=False)+'\n')
