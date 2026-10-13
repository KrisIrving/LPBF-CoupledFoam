"""M2A-02C native static joint wall/pressure candidate; unchanged14-case matrix."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

spec=importlib.util.spec_from_file_location('joint_surface_baseline',Path(__file__).with_name('prepare-m2a-surface.py'))
surface=importlib.util.module_from_spec(spec);spec.loader.exec_module(surface)
CASES=surface.CASES


def prepare(root,selection,ranks=None,spatial_only=False,reconstruction='trilinear',communication_audit=False,control_volume_audit=False):
    if reconstruction not in ('trilinear','fluidCubicAdaptive','completeWallBasis'):raise ValueError('Unknown reconstruction')
    if spatial_only and selection in ('fixed-restart','fixed-mpi2'):
        raise ValueError('Spatial reference uses twelve geometry cases only')
    if ranks is not None and ranks not in (1,2,4,8):
        raise ValueError('Joint ranks must be 1, 2, 4 or 8')
    if ranks is not None and selection=='fixed-restart' and ranks!=1:
        raise ValueError('Common restart remains serial')
    if ranks is not None and selection=='fixed-mpi2' and ranks!=2:
        raise ValueError('MPI comparison requires 2 ranks')
    meta=surface.prepare(root,selection)
    if ranks is not None:
        split={1:(1,1,1),2:(2,1,1),4:(2,2,1),8:(2,2,2)}[ranks]
        (root/'system/decomposeParDict').write_text(
            surface.baseline.foam_header('dictionary','decomposeParDict')+
            f'numberOfSubdomains {ranks};\nmethod simple;\n'
            f'simpleCoeffs {{ n ({split[0]} {split[1]} {split[2]}); delta 0.001; }}\n')
        meta['execution']={'ranks':ranks,'decomposition':list(split)}
        meta['nprocs']=ranks
    p=root/'constant/mechanicalProperties'
    text=p.read_text().split('constraintScheme surfaceExtension;')[0]
    p.write_text(text+'''constraintScheme jointSurface;
surfaceTargetTolerance 1e-7;
continuityTolerance 1e-7;
boundaryTreatment compatibleGauss;
jointKrylovBudget 64;
jointWarmStart true;
jointFullPredictor false;
jointOperatorAudit true;
jointOuterMix true;
jointOuterHistory 4;
jointSparseTransfer true;
''')
    solution=root/'system/fvSolution'
    old='p { solver PCG; preconditioner DIC; tolerance 1e-14; relTol 0; }'
    replacement='''p
    {
        solver GAMG;
        tolerance 1e-14;
        relTol 0;
        smoother DICGaussSeidel;
        cacheAgglomeration true;
        agglomerator faceAreaPair;
        nCellsInCoarsestLevel 50;
        nPreSweeps 0;
        nPostSweeps 2;
        nFinestSweeps 2;
        mergeLevels 1;
        maxIter 1000;
        jointPressurePolish true;
        jointGAMGIterations 80;
    }'''
    text=solution.read_text()
    if old not in text:raise ValueError('Expected baseline pressure controls not found')
    solution.write_text(text.replace(old,replacement))
    marker_count=math.ceil(4*math.pi*meta['radius']**2/(2.25*(.12/meta['mesh_n'])**2))
    meta['package']='M2A-02C'
    meta['surface_controls']={'scheme':'jointSurface','max_outer':1,'target_defect_m_s':1e-7}
    meta['joint_controls']={'markers':marker_count,'min_rank_pivot':1e-8,'krylov_budget':64,
        'force_exchange_N':1e-12,'torque_exchange_N_m':1e-14,'work_exchange_W':1e-14,
        'boundary_net_flux_m3_s':1e-18,'boundary_correction_m_s':1e-10}
    meta['joint_algorithm']={'revision':8,'warm_start':True,'full_predictor':False,
                             'response_reuse':True,'pressure_solver':'GAMG',
                             'pressure_tolerance':1e-14,'pressure_relTol':0.}
    meta['pressure_polish']={'gamg_iterations':80,'total_iterations':1000,
                            'tolerance':1e-14,'fallback':'PCG/DIC'}
    meta['outer_mix']={'enabled':True,'factor_max':8.,'history_depth':4,'kind':'residual_subspace',
                       'divergence_per_s':1e-7,
                       'accept':'actual momentum L1 decreases; original wall/div gates'}
    meta['transfer']={'sparse':True,'reference_relative':1e-12,'scaling':'once_per_frozen_matrix',
                      'profiling':'rank0 response phases; includes first-window response probes'}
    meta['operator_audit']={'version':1,'relative_tolerance':1e-8,
                            'divergence_per_s':1e-7,'manufactured_residual_m_s':2e-8,
                            'divergence_speed_factor':1e6,
                            'full_impulse':'diagnostic_only_diagonal_response',
                            'normal_gain':'diagnostic_only_not_full_rank_proof'}
    meta['scope']='Static single noncontact sphere; paired trilinear J and integrated-force transpose spread. '
    meta['scope']+='FV face-flux pressure-projected GMRES correction plus residual-controlled PISO. '
    meta['scope']+='No ghost penalty. q4 inventory retained; no moving GCL/contact/heat/scalable mapping validation.'
    if control_volume_audit:
        if spatial_only:raise ValueError('Control volume audit requires time advancement')
        with p.open('a') as stream:stream.write('controlVolumeAudit true;\n')
        meta['control_volume_audit']={'padding_h':[6,8],'matrix_match_impulse_N_s':1e-12,
            'face_match_impulse_N_s':1e-12,'angular_match_impulse_N_m_s':meta['radius']*1e-12,
            'source_force_N':1e-12,'source_torque_N_m':1e-14,
            'symmetric_stress_audit':True,
            'scope':'Static Euler/Gauss-linear mesh-face ledger; torque arm correction explicit; diagnostic only'}
    if communication_audit:
        with p.open('a') as stream:stream.write('jointCommunicationAudit true;\n')
        meta['communication_audit']={'enabled':True,'scope':'replicated response stage/input/cache agreement; per-rank last-stage breadcrumb',
                                     'native_validation':'pending','changes_numerical_operator':False}
    if spatial_only or reconstruction!='trilinear':
        reference_spec=importlib.util.spec_from_file_location('spatial_reference',Path(__file__).with_name('audit-m2a-spatial-consistency.py'))
        reference=importlib.util.module_from_spec(reference_spec);reference_spec.loader.exec_module(reference)
        expected=reference.fixture(meta['mode'],meta['mesh_n'],meta['centre_phase_h'])
        if reconstruction=='completeWallBasis':
            artifact=json.loads((Path(__file__).parents[1]/'docs/m2a-wall-basis-audit.json').read_text())
            geometry=next(c for c in artifact['cases'] if c['mesh_n']==meta['mesh_n'] and c['phase_h']==meta['centre_phase_h'])
            if not artifact['all_geometry_admitted']:raise ValueError('Incomplete wall basis geometry admission')
            expected.update(markers=geometry['modes'],wall_rms_m_s=geometry['references'][meta['mode']]['wall_rms_m_s'],
                            wall_max_m_s=geometry['references'][meta['mode']]['wall_max_m_s'])
            meta['joint_controls']['markers']=geometry['modes']
            meta['complete_wall_basis']={'geometry':geometry,'coverage_tolerance':1e-10,'orthogonality_tolerance':1e-10,'reference':expected}
            with p.open('a') as stream:
                stream.write(f'jointReconstruction completeWallBasis;\nexpectedWallModes {geometry["modes"]};\n')
                stream.write(f'referenceWallRms {expected["wall_rms_m_s"]:.17g};\nreferenceWallMax {expected["wall_max_m_s"]:.17g};\nreferenceStressRatio {expected["stress_load_ratio"]:.17g};\n')
            meta['scope']=meta['scope'].replace('paired trilinear J','orthonormal complete wall J')
            meta['scope']+=' Modal loads use basis zeroth/first moments, not point-force interpretation.'
        elif reconstruction!='trilinear':
            fluid_spec=importlib.util.spec_from_file_location('fluid_reference',Path(__file__).with_name('audit-m2a-fluid-stencil.py'))
            fluid=importlib.util.module_from_spec(fluid_spec);fluid_spec.loader.exec_module(fluid)
            geometry=fluid.fixture(meta['mesh_n'],meta['centre_phase_h'],3,True)
            expected['wall_rms_m_s']=geometry[meta['mode']+'_wall_rms_m_s']
            expected['wall_max_m_s']=geometry[meta['mode']+'_wall_max_m_s']
            meta['joint_reconstruction']={'scheme':reconstruction,'moment_tolerance':1e-10,'coefficient_l1_max':8,'geometry':geometry,'reference':expected}
            with p.open('a') as stream:stream.write('jointReconstruction fluidCubicAdaptive;\n')
            with p.open('a') as stream:
                stream.write(f'referenceWallRms {expected["wall_rms_m_s"]:.17g};\nreferenceWallMax {expected["wall_max_m_s"]:.17g};\nreferenceStressRatio {expected["stress_load_ratio"]:.17g};\n')
            meta['scope']=meta['scope'].replace('paired trilinear J','paired fluid-side cubic J')
        if spatial_only:meta['spatial_reference']=expected
    if spatial_only:
        surface.baseline.control(root,end=0)
        with p.open('a') as stream:stream.write('spatialAuditOnly true;\n')
        meta['scope']='Analytic-field native J/stress consistency only; no CFD/DEM time advancement or physical pass.'
    (root/'M2A_META.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta


def one_window(root):
    meta=json.loads((root/'M2A_META.json').read_text())
    if 'control_volume_audit' not in meta or meta.get('restart'):
        raise ValueError('One-window preflight requires a fresh control-volume case')
    surface.baseline.control(root,end=meta['delta_t'])
    meta['expected_rows']=1
    meta['preflight_only']=True
    (root/'M2A_META.json').write_text(json.dumps(meta,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',type=Path)
    parser.add_argument('--selection',choices=CASES)
    parser.add_argument('--ranks',type=int,choices=(1,2,4,8))
    parser.add_argument('--spatial-only',action='store_true')
    parser.add_argument('--reconstruction',choices=('trilinear','fluidCubicAdaptive','completeWallBasis'),default='trilinear')
    parser.add_argument('--communication-audit',action='store_true')
    parser.add_argument('--control-volume-audit',action='store_true')
    parser.add_argument('--one-window',action='store_true')
    parser.add_argument('--half',action='store_true')
    parser.add_argument('--restart',action='store_true')
    args=parser.parse_args()
    if args.one_window:one_window(args.case)
    elif args.restart:
        meta=json.loads((args.case/'M2A_META.json').read_text())
        surface.baseline.control(args.case,start=.002)
        surface.baseline.write_dem(args.case,meta['mode'],restart=True)
    elif args.half:surface.baseline.control(args.case,end=.002)
    elif args.selection:prepare(args.case,args.selection,args.ranks,args.spatial_only,args.reconstruction,args.communication_audit,args.control_volume_audit)
    else:parser.error('Provide selection, half or restart')
