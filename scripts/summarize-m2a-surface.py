"""Predeclared surface-candidate gates; retains every failed/incomplete case."""
import argparse
import importlib.util
import json
import math
import re
from pathlib import Path


def load(name,filename):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(filename))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


baseline=load('mechanical_audit','summarize-m2a.py')
prepare=load('surface_prepare','prepare-m2a-surface.py')


def evaluate(case):
    entry,rows,meta=baseline.evaluate(case)
    logs=('log.first','log.restart') if meta['restart'] else ('log.solver',)
    for filename in logs:
        text=(case/filename).read_text(errors='replace')
        scheme=meta['surface_controls'].get('scheme','surfaceExtension')
        if 'M2A_CONSTRAINT scheme='+scheme not in text:
            raise ValueError('Missing native '+scheme+' scheme marker')
        if scheme=='jointSurface':
            if 'M2A_JOINT markers=' not in text or 'M2A_BOUNDARY treatment=compatibleGauss' not in text:
                raise ValueError('Missing native joint/compatible boundary marker')
            match=re.search(r'M2A_BOUNDARY treatment=compatibleGauss rawFlux=(\S+) gaussFlux=(\S+) finalFlux=(\S+) normalCorrection=(\S+)',text)
            if not match or any(not math.isfinite(float(x)) for x in match.groups()):
                raise ValueError('Invalid native boundary compatibility diagnostics')
            if abs(float(match.group(3)))>meta['joint_controls']['boundary_net_flux_m3_s'] or abs(float(match.group(4)))>meta['joint_controls']['boundary_correction_m_s']:
                raise ValueError('Native boundary compatibility threshold exceeded')
            if 'joint_algorithm' in meta:
                marker=re.search(r'M2A_JOINT_ALGORITHM warmStart=(1|true) fullPredictor=(1|true) responseReuse=1',text)
                if not marker or 'GAMG:  Solving for p,' not in text:
                    raise ValueError('Missing accelerated joint algorithm/GAMG native evidence')
        if meta['surface_controls'].get('reconstruction')=='quadratic' and 'M2A_RECONSTRUCTION order=quadratic' not in text:
            raise ValueError('Missing native quadratic reconstruction marker')
    controls=meta['surface_controls']
    metrics,checks=entry['metrics'],entry['checks']
    metrics['surface_target_defect_m_s']=max(abs(r['surface_target_defect']) for r in rows)
    metrics['surface_correctors_max']=max(r['surface_correctors'] for r in rows)
    metrics['surface_correctors_mean']=sum(r['surface_correctors'] for r in rows)/len(rows)
    metrics['pressure_correctors_total']=sum(r['pressure_correctors_total'] for r in rows)
    metrics['window_wall_seconds_total']=sum(r['window_wall_seconds'] for r in rows)
    checks['surface_target_convergence']=metrics['surface_target_defect_m_s']<=controls['target_defect_m_s']
    checks['surface_iteration_budget']=all(r['surface_correctors']==int(r['surface_correctors'])
        and 1<=r['surface_correctors']<=controls['max_outer'] for r in rows)
    checks['total_pressure_counts']=all(r['pressure_correctors_total']==int(r['pressure_correctors_total'])
        and r['surface_correctors']*meta['solver_controls']['min_correctors']
        <=r['pressure_correctors_total']<=r['surface_correctors']*meta['solver_controls']['max_correctors'] for r in rows)
    checks['positive_window_cost']=all(r['window_wall_seconds']>0 for r in rows)
    if 'continuity_per_s' in meta['solver_controls']:
        checks['continuity_convergence']=metrics['divergence']<=meta['solver_controls']['continuity_per_s']
    if 'joint_controls' in meta:
        joint=meta['joint_controls']
        metrics['joint_wall_residual']=max(abs(r['joint_wall_residual']) for r in rows)
        metrics['joint_iterations_total']=sum(r['joint_iterations'] for r in rows)
        metrics['joint_pressure_solves_total']=sum(r['joint_pressure_solves'] for r in rows)
        checks['joint_wall_constraint']=metrics['joint_wall_residual']<=controls['target_defect_m_s']
        checks['joint_marker_rank']=all(r['joint_markers']==joint['markers'] and r['joint_min_pivot']>joint['min_rank_pivot'] for r in rows)
        checks['joint_iteration_counts']=all(r['joint_iterations']==int(r['joint_iterations']) and 0<=r['joint_iterations']<=r['pressure_correctors']*joint['krylov_budget']
            and r['joint_pressure_solves']==int(r['joint_pressure_solves']) and 0<=r['joint_pressure_solves']<=r['pressure_correctors']*(joint['krylov_budget']+2) for r in rows)
        for name,limit in (('force','force_exchange_N'),('torque','torque_exchange_N_m'),('work','work_exchange_W')):
            metrics['joint_'+name+'_exchange_error']=max(abs(r['joint_'+name+'_exchange_error']) for r in rows)
            checks['joint_'+name+'_exchange']=metrics['joint_'+name+'_exchange_error']<=joint[limit]
        if meta['restart'] and meta.get('joint_algorithm',{}).get('warm_start'):
            checkpoint=(case/'0.002/couplingState').read_text()
            match=re.search(r'jointLoadSeed\s+(?:(\d+)\s*)?\(([^)]*)\)\s*;',checkpoint,re.S)
            uniform=re.search(r'jointLoadSeed\s+(\d+)\s*\{\s*([^}]+)\}\s*;',checkpoint,re.S)
            size=3*joint['markers']
            if match:
                seed=[float(x) for x in match.group(2).split()]
                declared=int(match.group(1)) if match.group(1) else len(seed)
            elif uniform:
                declared=int(uniform.group(1))
                seed=[float(uniform.group(2))] if declared==size else []
            else:raise ValueError('Missing joint warm-start load seed')
            checks['joint_warm_checkpoint']=bool(
                declared==size and (len(seed)==size if match else len(seed)==1)
                and all(math.isfinite(x) for x in seed)
                and re.search(r'jointWarmStart\s+(?:1|true|on)\s*;',checkpoint)
                and re.search(r'jointFullPredictor\s+(?:1|true|on)\s*;',checkpoint))
    radius=meta['radius']
    reference=(6*math.pi*meta['rho']*meta['nu']*radius*meta['reference_flow_speed']
        if meta['mode']=='fixed' else 8*math.pi*meta['rho']*meta['nu']*radius**3*meta['reference_rotation_speed'])
    prefix='f' if meta['mode']=='fixed' else 't'
    stress=[sum(r['stress_pressure_'+prefix+c]+r['stress_viscous_'+prefix+c]
                for r in rows[-10:])/10 for c in 'xyz']
    ledger=[sum(r[prefix+c] for r in rows[-10:])/10 for c in 'xyz']
    metrics['stress_ledger_relative']=baseline.norm([a-b for a,b in zip(stress,ledger)])/reference
    metrics['stress_stokes_relative']=abs((stress[0] if prefix=='f' else -stress[2])/reference-1)
    metrics['mean_reference_load_ratio']=(ledger[0] if prefix=='f' else -ledger[2])/reference
    if meta['mesh_n']==64:
        checks['fine64_reference']=metrics['stokes_relative_error']<=meta['surface_thresholds']['fine64_reference_relative']
        checks['fine64_stress_ledger']=metrics['stress_ledger_relative']<=meta['surface_thresholds']['fine64_stress_ledger_relative']
    entry['passed']=all(checks.values())
    return entry,rows,meta


def summarize(root,package='M2A-02B'):
    result={'package':package,'cases':{},'comparisons':{},
            'scope':'Static noncontact single-sphere surface candidate only. '
                    'No moving GCL, general angular conservation, contact, heat or LPBF validation.'}
    histories,metadata={},{}
    for name in prepare.CASES:
        try:
            result['cases'][name],histories[name],metadata[name]=evaluate(root/name)
        except (OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as error:
            result['cases'][name]={'passed':False,'error':str(error),
                                  'diagnostics':baseline.failure_diagnostics(root/name)}
    for name in ('fixed-restart','fixed-mpi2'):
        if 'fixed-coarse' in histories and name in histories:
            result['comparisons'][name]=baseline.compare(histories['fixed-coarse'],histories[name],
                                                         metadata['fixed-coarse']['thresholds'])
        else:
            result['comparisons'][name]={'passed':False,'status':'not_evaluated'}
    for mode in ('fixed','rotate'):
        for offset in ('','-offset'):
            for coarse,fine in (('coarse','fine'),('fine','finer')):
                a,b=mode+'-'+coarse+offset,mode+'-'+fine+offset
                label=f'{mode}-{coarse}-to-{fine}{offset}'
                am,bm=result['cases'][a].get('metrics',{}),result['cases'][b].get('metrics',{})
                if 'stokes_relative_error' in am and 'stokes_relative_error' in bm:
                    increase=bm['stokes_relative_error']-am['stokes_relative_error']
                    result['comparisons'][label]={'passed':increase<=metadata[b]['thresholds']['refinement_error_increase'],
                        'error_increase':increase,'interpretation':'No asymptotic order claim'}
                else:
                    result['comparisons'][label]={'passed':False,'status':'not_evaluated'}
        for level in ('coarse','fine','finer'):
            a,b=mode+'-'+level,mode+'-'+level+'-offset'
            am,bm=result['cases'][a].get('metrics',{}),result['cases'][b].get('metrics',{})
            label=mode+'-'+level+'-phase'
            if 'mean_reference_load_ratio' in am and 'mean_reference_load_ratio' in bm:
                difference=abs(am['mean_reference_load_ratio']-bm['mean_reference_load_ratio'])
                result['comparisons'][label]={'passed':difference<=metadata[b]['surface_thresholds']['phase_load_difference_relative'],
                                              'load_difference_relative':difference}
            else:
                result['comparisons'][label]={'passed':False,'status':'not_evaluated'}
    result['passed']=all(x['passed'] for section in ('cases','comparisons') for x in result[section].values())
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=summarize(args.report)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for section in ('cases','comparisons'):
        for name,entry in result[section].items():
            print(name,'PASS' if entry['passed'] else 'FAIL',entry.get('error',''))
    raise SystemExit(0 if result['passed'] else 1)
