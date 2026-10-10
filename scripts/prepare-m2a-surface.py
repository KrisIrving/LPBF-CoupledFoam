"""M2A-02B static surface candidate: three grids, two phases, MPI/restart."""
import argparse
import importlib.util
import json
from pathlib import Path

spec = importlib.util.spec_from_file_location('prepare_baseline',Path(__file__).with_name('prepare-m2a.py'))
baseline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(baseline)
CASES = {f'{mode}-{level}{offset}':(mode,n,phase)
         for mode in ('fixed','rotate')
         for level,n in (('coarse',36),('fine',48),('finer',64))
         for offset,phase in (('',0.),('-offset',.25))}
CASES.update({'fixed-restart':('fixed',36,0.),'fixed-mpi2':('fixed',36,0.)})


def prepare(root,selection):
    mode,n,phase = CASES[selection]
    # Reuse baseline FV/DEM controls, not a manufactured particle force.
    meta = baseline.prepare(root,mode+'-coarse')
    h=.12/n
    mesh=root/'system/blockMeshDict'
    mesh.write_text(mesh.read_text().replace('(36 36 36)',f'({n} {n} {n})'))
    centre=[phase*h]*3
    data=root/'particle.data'
    data.write_text(data.read_text().replace('7 1 0.02 1000 0.0 0 0',
                    '7 1 0.02 1000 '+' '.join(format(x,'.17g') for x in centre)))
    properties=root/'constant/mechanicalProperties'
    properties.write_text(properties.read_text()+'''
constraintScheme surfaceExtension;
maxSurfaceCorrectors 32;
surfaceTargetTolerance 1e-7;
''')
    meta.update(package='M2A-02B',selection=selection,mesh_n=n,
                initial_centre=centre,centre_phase_h=phase,
                nprocs=2 if selection=='fixed-mpi2' else 1,
                restart=selection=='fixed-restart')
    meta['surface_controls']={'scheme':'surfaceExtension','max_outer':32,
                              'target_defect_m_s':1e-7,'stress_z_nodes':12,'stress_phi_nodes':24}
    meta['surface_thresholds']={'fine64_reference_relative':.10,
                                'fine64_stress_ledger_relative':.10,
                                'phase_load_difference_relative':.02}
    meta['scope']='Static single noncontact sphere. Cell-centre interior constraint with true-radius '
    meta['scope']+='fluid-side extension; q4 fictitious-fluid volume retained. Wall interpolation slip '
    meta['scope']+='and one-sided stress diagnostic. No moving geometry/contact/heat/GCL validation.'
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
        baseline.control(args.case,start=.002)
        baseline.write_dem(args.case,meta['mode'],restart=True)
    elif args.half:
        baseline.control(args.case,end=.002)
    elif args.selection:
        prepare(args.case,args.selection)
    else:
        parser.error('Provide selection, half or restart')
