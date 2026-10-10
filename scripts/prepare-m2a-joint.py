"""M2A-02C native static joint wall/pressure candidate; unchanged14-case matrix."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

spec=importlib.util.spec_from_file_location('joint_surface_baseline',Path(__file__).with_name('prepare-m2a-surface.py'))
surface=importlib.util.module_from_spec(spec);spec.loader.exec_module(surface)
CASES=surface.CASES


def prepare(root,selection):
    meta=surface.prepare(root,selection)
    p=root/'constant/mechanicalProperties'
    text=p.read_text().split('constraintScheme surfaceExtension;')[0]
    p.write_text(text+'''constraintScheme jointSurface;
surfaceTargetTolerance 1e-7;
continuityTolerance 1e-7;
boundaryTreatment compatibleGauss;
jointKrylovBudget 64;
''')
    marker_count=math.ceil(4*math.pi*meta['radius']**2/(2.25*(.12/meta['mesh_n'])**2))
    meta['package']='M2A-02C'
    meta['surface_controls']={'scheme':'jointSurface','max_outer':1,'target_defect_m_s':1e-7}
    meta['joint_controls']={'markers':marker_count,'min_rank_pivot':1e-8,'krylov_budget':64,
        'force_exchange_N':1e-12,'torque_exchange_N_m':1e-14,'work_exchange_W':1e-14,
        'boundary_net_flux_m3_s':1e-18,'boundary_correction_m_s':1e-10}
    meta['scope']='Static single noncontact sphere; paired trilinear J and integrated-force transpose spread. '
    meta['scope']+='FV face-flux pressure-projected GMRES correction plus residual-controlled PISO. '
    meta['scope']+='No ghost penalty. q4 inventory retained; no moving GCL/contact/heat/scalable mapping validation.'
    (root/'M2A_META.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',type=Path)
    parser.add_argument('--selection',choices=CASES)
    parser.add_argument('--half',action='store_true')
    parser.add_argument('--restart',action='store_true')
    args=parser.parse_args()
    if args.restart:
        meta=json.loads((args.case/'M2A_META.json').read_text())
        surface.baseline.control(args.case,start=.002)
        surface.baseline.write_dem(args.case,meta['mode'],restart=True)
    elif args.half:surface.baseline.control(args.case,end=.002)
    elif args.selection:prepare(args.case,args.selection)
    else:parser.error('Provide selection, half or restart')
