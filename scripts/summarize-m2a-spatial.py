"""Compare native analytic-field diagnostics with independent Python oracle."""
import argparse
import json
import math
import re
from pathlib import Path


def reference_values(text,expected):
    matches=re.findall(r'M2A_SPATIAL_REFERENCE markers=(\d+) wallRms=(\S+) wallMax=(\S+) stressRatio=(\S+)',text)
    if len(matches)!=1:raise ValueError('Missing or duplicate native spatial reference')
    count,rms,maximum,stress=matches[0];actual={'wall_rms_m_s':float(rms),'wall_max_m_s':float(maximum),'stress_load_ratio':float(stress)}
    if not all(math.isfinite(value) for value in actual.values()):
        raise ValueError('Nonfinite native spatial reference')
    checks={'markers':int(count)==expected['markers']}
    for key,value in actual.items():
        scale=expected['reference_speed_m_s'] if key.startswith('wall_') else 1.
        checks[key]=math.isfinite(value) and abs(value-expected[key])/scale<=1e-8
    return {'passed':all(checks.values()),'checks':checks,'actual':actual,'expected':expected}


def reconstruction_values(text,controls):
    matches=re.findall(r'M2A_RECONSTRUCTION scheme=(\w+) moment=(\S+) l1=(\S+) supportMin=(\S+) supportMax=(\S+) expanded=(\S+)',text)
    if len(matches)!=1 or matches[0][0]!=controls['scheme']:raise ValueError('Missing reconstruction evidence')
    moment,l1,minimum,maximum,expanded=map(float,matches[0][1:]);geometry=controls['geometry']
    if (not all(math.isfinite(x) for x in (moment,l1,minimum,maximum,expanded))
        or not 0<=moment<=controls['moment_tolerance'] or not 0<l1<=controls['coefficient_l1_max']
        or (minimum,maximum,expanded)!=(geometry['support_min'],geometry['support_max'],geometry['expanded_markers'])):
        raise ValueError('Failed reconstruction geometry/moment/coefficient evidence')
    return {'moment':moment,'l1':l1,'support_min':minimum,'support_max':maximum,'expanded':expanded}


def basis_values(text,controls):
    matches=re.findall(r'M2A_WALL_BASIS nodes=(\S+) modes=(\S+) coverage=(\S+) orthogonality=(\S+) signature=(b\w+)',text)
    if len(matches)!=1:raise ValueError('Missing complete wall basis evidence')
    nodes,modes,coverage,orthogonality=map(float,matches[0][:4]);geometry=controls['geometry']
    if (not all(math.isfinite(x) for x in (nodes,modes,coverage,orthogonality))
        or nodes!=geometry['full_wall_nodes'] or modes!=geometry['modes']
        or not 0<=coverage<=controls['coverage_tolerance']
        or not 0<=orthogonality<=controls['orthogonality_tolerance']):
        raise ValueError('Complete wall basis coverage/rank/orthogonality failed')
    return {'nodes':nodes,'modes':modes,'coverage':coverage,'orthogonality':orthogonality,'signature':matches[0][4]}


def evaluate(case):
    meta=json.loads((case/'M2A_META.json').read_text());expected=meta['spatial_reference']
    text=(case/'log.solver').read_text(errors='replace')
    if (case/'result.txt').read_text().strip()!='exit_status=0' or 'M2A_EXECUTION_COMPLETE' not in text:
        raise ValueError('Native reference did not complete')
    rank=re.search(r'^nProcs\s*:\s*(\d+)\s*$',text,re.M)
    if not rank or int(rank[1])!=meta['execution']['ranks']:raise ValueError('Native MPI ranks mismatch')
    result=reference_values(text,expected)
    if 'joint_reconstruction' in meta:result['reconstruction']=reconstruction_values(text,meta['joint_reconstruction'])
    if 'complete_wall_basis' in meta:result['complete_wall_basis']=basis_values(text,meta['complete_wall_basis'])
    return result


def summarize(root):
    cases={}
    for mode in ('fixed','rotate'):
        for level in ('coarse','fine','finer'):
            for offset in ('','-offset'):
                name=f'{mode}-{level}{offset}'
                try:cases[name]=evaluate(root/name)
                except (OSError,KeyError,ValueError) as error:cases[name]={'passed':False,'error':str(error)}
    return {'scope':__doc__+' No CFD convergence, load accuracy or M2A closure.',
            'cases':cases,'passed':all(c['passed'] for c in cases.values())}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('report',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=summarize(args.report)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for name,c in result['cases'].items():print(name,'REFERENCE_MATCH' if c['passed'] else 'FAIL',c.get('error',''))
    raise SystemExit(0 if result['passed'] else 1)
