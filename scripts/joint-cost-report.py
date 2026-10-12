"""Diagnostic progress/cost of completed windows, including failed runs; no gates."""
import argparse
import csv
import json
import math
import re
from pathlib import Path


def collect(root):
    result={'scope':'Diagnostic only; incomplete trajectories cannot pass physical gates. '
                     'Rank0 timings exclude unfinished windows and other work.', 'cases':{}}
    for case in sorted(root.iterdir()):
        if not case.is_dir() or not (case/'M2A_META.json').exists():continue
        rows=[]
        history=case/'mechanical-history.csv'
        if history.exists():
            with history.open() as stream:rows=list(csv.DictReader(stream))
        entry={'completed_windows':len(rows),
               'execution':(case/'result.txt').read_text().strip() if (case/'result.txt').exists() else 'missing',
               'logs':{}}
        for filename in ('log.solver','log.first','log.restart'):
            path=case/filename
            if not path.exists():continue
            text=path.read_text(errors='replace');costs=[];invalid=0
            keys=('response','spread','projection','interpolation','scaling')
            for line in re.findall(r'M2A_RESPONSE_COST [^\n]+',text):
                try:
                    values=dict(re.findall(r'(\w+)=(\S+)',line))
                    cost={key:float(values[key]) for key in keys}
                    if not all(math.isfinite(v) and v>=0 for v in cost.values()):raise ValueError()
                    costs.append(cost)
                except (KeyError,ValueError):invalid+=1
            totals={key:sum(c[key] for c in costs) for key in keys}
            ranks=re.search(r'^nProcs\s*:\s*(\d+)\s*$',text,re.M)
            entry['logs'][filename]={'native_ranks':int(ranks[1]) if ranks else None,
                'cost_windows':len(costs),'invalid_cost_records':invalid,'cost_seconds':totals,
                'projection_fraction_of_response':totals['projection']/totals['response'] if totals['response'] else None,
                'explicit_failures':re.findall(r'M2A_FAIL[^\n]+',text)}
        result['cases'][case.name]=entry
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(collect(args.report),indent=2,allow_nan=False)+'\n')
