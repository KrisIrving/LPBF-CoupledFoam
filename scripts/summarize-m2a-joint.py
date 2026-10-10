"""Same physical/refinement/phase gates plus joint operator checks."""
import argparse
import importlib.util
import json
from pathlib import Path

spec=importlib.util.spec_from_file_location('joint_common_gates',Path(__file__).with_name('summarize-m2a-surface.py'))
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


def summarize(root):
    result=audit.summarize(root,package='M2A-02C')
    # Do not let a successful legacy surface run masquerade as a joint run.
    for name,entry in result['cases'].items():
        try:
            meta=json.loads((root/name/'M2A_META.json').read_text())
            if meta['package']!='M2A-02C' or meta['surface_controls']['scheme']!='jointSurface' or 'joint_controls' not in meta:
                raise ValueError('Wrong package/scheme for joint verification')
        except (OSError,ValueError,KeyError) as error:
            entry['passed']=False;entry['error']=str(error)
    result['passed']=all(x['passed'] for section in ('cases','comparisons') for x in result[section].values())
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=summarize(args.report)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for section in ('cases','comparisons'):
        for name,entry in result[section].items():print(name,'PASS' if entry['passed'] else 'FAIL',entry.get('error',''))
    raise SystemExit(0 if result['passed'] else 1)
