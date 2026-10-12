"""Generator/gate regressions using synthetic histories, not CFD evidence."""
import importlib.util
import json
import re
from pathlib import Path
import tempfile
import unittest

root=Path(__file__).parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

prepare=load('joint_package_prepare',root/'scripts/prepare-m2a-joint.py')
summary=load('joint_package_summary',root/'scripts/summarize-m2a-joint.py')
geometry=load('joint_package_geometry',root/'scripts/audit-m2a-joint-stencils.py')
cost_report=load('joint_cost_report',root/'scripts/joint-cost-report.py')
synthetic=load('joint_legacy_fixture',root/'tests/test_m2a_surface.py')


class JointPackageTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.writer=synthetic.SurfaceGateTests();self.writer.root=self.root
    def tearDown(self):self.temp.cleanup()
    def fixture(self,name):
        case,rows,_=self.writer.fixture(name)
        meta=prepare.prepare(case,name)
        if meta['restart']:
            size=3*meta['joint_controls']['markers']
            (case/'0.002/couplingState').write_text('jointWarmStart true;\njointFullPredictor false;\njointOuterMix true;\njointOuterHistory 4;\njointSparseTransfer true;\njointLoadSeed '
                +str(size)+'\n('+ ' '.join('0' for _ in range(size))+');\n')
        for log in ('log.first','log.restart') if meta['restart'] else ('log.solver',):
            (case/log).write_text('M2A_CONSTRAINT scheme=jointSurface\nM2A_JOINT markers='+str(meta['joint_controls']['markers'])+' minPivot=.9\n'
                'M2A_BOUNDARY treatment=compatibleGauss rawFlux=-4e-14 gaussFlux=7e-18 finalFlux=1e-22 normalCorrection=1e-16\n'
                'M2A_JOINT_ALGORITHM warmStart=1 fullPredictor=0 responseReuse=1\nGAMG:  Solving for p,\n'
                'M2A_PRESSURE_SOLVE gamg=80 pcg=12 residual=9e-15 tolerance=1e-14 budget=1000\n'
                'M2A_OUTER_MODE enabled=1\n'
                'M2A_OUTER_HISTORY depth=4\n'
                'M2A_TRANSFER sparse=1\nM2A_TRANSFER_AUDIT interpolation=1e-15 spreading=1e-15 scaling=1e-15\n'
                'M2A_RESPONSE_COST time=.0001 response=1 spread=.1 projection=.8 interpolation=.05 scaling=.01\n'
                'M2A_OUTER_MIX factor=2 accepted=1 before=2e-10 after=1e-10 wall=1e-8 div=1e-12 linear=2e-12 history=3\n'
                'M2A_OPERATOR_AUDIT time=.0001 zero=0 repeat=0 linear=1e-12 scale=1e-12 div=1e-12 gauge=0 diagonal=1e-15 fullImpulse=1e-6 normalGain=.1 manufactured=1e-9 iterations=12 pressureSolves=21 seconds=.1\n'
                'M2A_OPERATOR_COMMIT time=.0001 source=1e-15 velocity=1e-15 pressure=1e-15 flux=1e-15 acceleration=1e-15\nM2A_EXECUTION_COMPLETE\n')
        for row in rows:
            row.update(surface_correctors=1,pressure_correctors_total=8,joint_iterations=8,
                joint_pressure_solves=24,joint_wall_residual=1e-8,joint_markers=meta['joint_controls']['markers'],
                joint_min_pivot=.9,joint_force_exchange_error=0.,joint_torque_exchange_error=0.,joint_work_exchange_error=0.)
        self.writer.write(case,rows);return case,rows,meta

    def test_full_matrix_generation_and_synthetic_gates(self):
        for name in prepare.CASES:
            case,rows,meta=self.fixture(name)
            properties=(case/'constant/mechanicalProperties').read_text()
            self.assertIn('constraintScheme jointSurface;',properties)
            self.assertIn('boundaryTreatment compatibleGauss;',properties)
            self.assertIn('jointWarmStart true;',properties)
            self.assertIn('jointFullPredictor false;',properties)
            self.assertIn('jointOperatorAudit true;',properties)
            self.assertIn('solver GAMG;',(case/'system/fvSolution').read_text())
            self.assertIn('tolerance 1e-14;',(case/'system/fvSolution').read_text())
            self.assertNotIn('constraintScheme surfaceExtension;',properties)
            self.assertEqual(meta['package'],'M2A-02C')
        self.assertTrue(summary.summarize(self.root)['passed'])

    def test_mpi_decomposition_and_native_rank_mismatch(self):
        case,rows,meta=self.fixture('fixed-finer')
        for ranks,split in ((1,'1 1 1'),(2,'2 1 1'),(4,'2 2 1'),(8,'2 2 2')):
            prepare.prepare(case,'fixed-finer',ranks)
            for row in rows:row['covered_ranks']=ranks
            self.writer.write(case,rows)
            self.assertIn('n ('+split+');',(case/'system/decomposeParDict').read_text())
            log=case/'log.solver';text=re.sub(r'^nProcs.*\n','',log.read_text(),flags=re.M)
            log.write_text('nProcs : '+str(ranks)+'\n'+text)
            self.assertTrue(summary.audit.evaluate(case)[0]['checks'].get('multi_rank_coverage',True))
            log.write_text('nProcs : 16\n'+text)
            with self.assertRaisesRegex(ValueError,'MPI rank evidence'):summary.audit.evaluate(case)
        with self.assertRaises(ValueError):prepare.prepare(case,'fixed-finer',3)
        with self.assertRaises(ValueError):prepare.prepare(case,'fixed-restart',8)
        with self.assertRaises(ValueError):prepare.prepare(case,'fixed-mpi2',8)

    def test_timeout_progress_is_diagnostic_and_nonfinite_cost_rejected(self):
        case,rows,meta=self.fixture('fixed-coarse')
        (case/'result.txt').write_text('exit_status=124\n')
        log=case/'log.solver';log.write_text(log.read_text()+'M2A_RESPONSE_COST response=nan spread=0 projection=0 interpolation=0 scaling=0\n')
        entry=cost_report.collect(self.root)['cases']['fixed-coarse']
        self.assertEqual(entry['completed_windows'],len(rows))
        self.assertEqual(entry['logs']['log.solver']['cost_windows'],1)
        self.assertEqual(entry['logs']['log.solver']['invalid_cost_records'],1)
        self.assertAlmostEqual(entry['logs']['log.solver']['projection_fraction_of_response'],.8)
        self.assertFalse(summary.summarize(self.root)['cases']['fixed-coarse']['passed'])

    def test_legacy_runs_cannot_masquerade_as_joint(self):
        for name in prepare.CASES:self.writer.fixture(name)
        result=summary.summarize(self.root)
        self.assertFalse(result['passed'])
        self.assertTrue(all(not entry['passed'] for entry in result['cases'].values()))

    def test_pressure_polish_rejects_missing_unconverged_and_excess_budget(self):
        case,_,_=self.fixture('fixed-coarse');p=case/'log.solver';original=p.read_text()
        for before,after in [('M2A_PRESSURE_SOLVE','OLD_PRESSURE_SOLVE'),
                             ('residual=9e-15','residual=1.5e-14'),
                             ('pcg=12','pcg=921'),('gamg=80','gamg=81'),
                             ('tolerance=1e-14','tolerance=2e-14')]:
            p.write_text(original.replace(before,after))
            with self.assertRaisesRegex(ValueError,'[Pp]ressure polish'):summary.audit.evaluate(case)

    def test_outer_mix_rejects_inconsistent_accepted_state_and_missing_mode(self):
        case,_,_=self.fixture('fixed-coarse');p=case/'log.solver';original=p.read_text()
        for before,after in [('M2A_OUTER_MODE enabled=1','M2A_OUTER_MODE enabled=0'),
                             ('factor=2','factor=9'),('after=1e-10','after=3e-10'),
                             ('linear=2e-12','linear=1e-5'),('wall=1e-8','wall=1e-5'),
                             ('div=1e-12','div=1e-3')]:
            p.write_text(original.replace(before,after))
            with self.assertRaisesRegex(ValueError,'[Oo]uter mixing'):summary.audit.evaluate(case)
        p.write_text(original.replace('accepted=1','accepted=0').replace('after=1e-10','after=3e-10'))
        self.assertTrue(summary.audit.evaluate(case)[0]['passed'])

    def test_outer_mix_checkpoint_required(self):
        case,_,_=self.fixture('fixed-restart');p=case/'0.002/couplingState'
        p.write_text(p.read_text().replace('jointOuterMix true;','jointOuterMix false;'))
        self.assertFalse(summary.audit.evaluate(case)[0]['passed'])

    def test_outer_history_depth_and_checkpoint_must_match(self):
        case,_,_=self.fixture('fixed-restart');p=case/'log.restart';original=p.read_text()
        for before,after in [('depth=4','depth=1'),('history=3','history=5')]:
            p.write_text(original.replace(before,after))
            with self.assertRaisesRegex(ValueError,'outer mixing'):summary.audit.evaluate(case)
        p.write_text(original);checkpoint=case/'0.002/couplingState'
        checkpoint.write_text(checkpoint.read_text().replace('jointOuterHistory 4;','jointOuterHistory 1;'))
        self.assertFalse(summary.audit.evaluate(case)[0]['passed'])

    def test_sparse_reference_and_cost_evidence_are_required(self):
        case,_,_=self.fixture('fixed-coarse');p=case/'log.solver';original=p.read_text()
        for before,after in [('M2A_TRANSFER sparse=1','M2A_TRANSFER sparse=0'),
                             ('M2A_TRANSFER_AUDIT','OLD_TRANSFER_AUDIT'),
                             ('spreading=1e-15','spreading=1e-5'),
                             ('M2A_RESPONSE_COST','OLD_RESPONSE_COST'),
                             ('projection=.8','projection=2'),('scaling=.01','scaling=nan')]:
            p.write_text(original.replace(before,after))
            with self.assertRaisesRegex(ValueError,'sparse transfer'):summary.audit.evaluate(case)

    def test_sparse_transfer_restart_mode_is_verified(self):
        case,_,_=self.fixture('fixed-restart');p=case/'0.002/couplingState'
        p.write_text(p.read_text().replace('jointSparseTransfer true;','jointSparseTransfer false;'))
        self.assertFalse(summary.audit.evaluate(case)[0]['passed'])

    def test_joint_failure_gates(self):
        case,rows,meta=self.fixture('fixed-coarse')
        rows[0].update(joint_wall_residual=2e-7,joint_min_pivot=1e-9,
            joint_iterations=513,joint_force_exchange_error=2e-12,
            joint_torque_exchange_error=2e-14,joint_work_exchange_error=2e-14)
        self.writer.write(case,rows)
        entry,_,_=summary.audit.evaluate(case)
        for name in ('joint_wall_constraint','joint_marker_rank','joint_iteration_counts',
                     'joint_force_exchange','joint_torque_exchange','joint_work_exchange'):
            self.assertFalse(entry['checks'][name])

    def test_native_boundary_and_joint_markers_required(self):
        case,rows,meta=self.fixture('fixed-coarse')
        p=case/'log.solver';text=p.read_text()
        p.write_text(text.replace('finalFlux=1e-22','finalFlux=1e-6'))
        with self.assertRaisesRegex(ValueError,'compatibility threshold'):summary.audit.evaluate(case)
        p.write_text(text.replace('M2A_JOINT markers=','M2A_OLD markers='))
        with self.assertRaisesRegex(ValueError,'joint/compatible'):summary.audit.evaluate(case)

    def test_six_geometric_stencils(self):
        result=geometry.run();self.assertTrue(result['passed'])
        self.assertEqual([case['markers'] for case in result['cases']],[51,51,90,90,159,159])
        self.assertGreater(min(case['min_normalised_cholesky_pivot'] for case in result['cases']),.9)

    def test_algorithm_evidence_cannot_be_skipped(self):
        case,rows,meta=self.fixture('fixed-coarse')
        p=case/'log.solver';p.write_text(p.read_text().replace('fullPredictor=0','fullPredictor=1'))
        with self.assertRaisesRegex(ValueError,'accelerated joint'):summary.audit.evaluate(case)

    def test_historical_revision2_mode_is_checked_against_its_metadata(self):
        case,rows,meta=self.fixture('fixed-restart')
        meta['joint_algorithm'].update(revision=2,full_predictor=True)
        meta.pop('operator_audit')
        (case/'M2A_META.json').write_text(json.dumps(meta))
        for name in ('log.first','log.restart'):
            p=case/name;p.write_text(p.read_text().replace('fullPredictor=0','fullPredictor=1'))
        p=case/'0.002/couplingState'
        p.write_text(p.read_text().replace('jointFullPredictor false;','jointFullPredictor true;'))
        self.assertTrue(summary.audit.evaluate(case)[0]['passed'])

    def test_operator_audit_evidence_and_failures(self):
        case,rows,meta=self.fixture('fixed-coarse');p=case/'log.solver';original=p.read_text()
        for before,after in [('M2A_OPERATOR_AUDIT','M2A_OLD_AUDIT'),('linear=1e-12','linear=1e-5'),
                             ('normalGain=.1','normalGain=nan'),('manufactured=1e-9','manufactured=1e-5'),
                             ('pressureSolves=21','pressureSolves=0'),('source=1e-15','source=1e-3'),
                             ('M2A_OPERATOR_COMMIT','M2A_OLD_COMMIT')]:
            p.write_text(original.replace(before,after))
            with self.assertRaisesRegex(ValueError,'joint operator audit'):summary.audit.evaluate(case)
        p.write_text(original.replace('fullImpulse=1e-6','fullImpulse=1'))
        self.assertTrue(summary.audit.evaluate(case)[0]['passed'])

    def test_restart_requires_its_own_operator_audit(self):
        case,rows,meta=self.fixture('fixed-restart');p=case/'log.restart'
        p.write_text(p.read_text().replace('M2A_OPERATOR_AUDIT','M2A_OLD_AUDIT'))
        with self.assertRaisesRegex(ValueError,'Missing actual joint operator audit'):summary.audit.evaluate(case)

    def test_failed_operator_evidence_is_retained_without_a_pass(self):
        case,rows,meta=self.fixture('fixed-coarse')
        (case/'result.txt').write_text('exit_status=1\n')
        with (case/'log.solver').open('a') as log:
            log.write('M2A_OPERATOR_PROBE time=.0001 linear=1e-3\n'
                      'M2A_FAIL rank=0: Actual joint operator probe failed before manufactured solve\n')
        entry=summary.summarize(self.root)['cases']['fixed-coarse']
        self.assertFalse(entry['passed'])
        self.assertIn('linear=1e-3',entry['diagnostics']['operator_logs']['log.solver']['probe'])
        self.assertTrue(entry['diagnostics']['operator_logs']['log.solver']['failures'])

    def test_warm_restart_seed_checked_independently(self):
        case,rows,meta=self.fixture('fixed-restart')
        checkpoint=case/'0.002/couplingState'
        text=checkpoint.read_text()
        checkpoint.write_text(text.replace('(0 ','(nan ',1))
        self.assertFalse(summary.audit.evaluate(case)[0]['checks']['joint_warm_checkpoint'])

        checkpoint.write_text(text.replace('jointWarmStart true;','jointWarmStart false;'))
        self.assertFalse(summary.audit.evaluate(case)[0]['checks']['joint_warm_checkpoint'])
        checkpoint.write_text('jointWarmStart true; jointFullPredictor false; jointLoadSeed 153{0};')
        self.assertTrue(summary.audit.evaluate(case)[0]['checks']['joint_warm_checkpoint'])
        checkpoint.write_text('jointWarmStart true; jointFullPredictor true; jointLoadSeed 153{0};')
        self.assertFalse(summary.audit.evaluate(case)[0]['checks']['joint_warm_checkpoint'])

    def test_adaptive_reconstruction_reference_mode_and_restart_evidence(self):
        for name in ('fixed-coarse','fixed-restart'):
            case,rows,_=self.fixture(name)
            meta=prepare.prepare(case,name,reconstruction='fluidCubicAdaptive')
            controls=meta['joint_reconstruction'];g=controls['geometry'];r=controls['reference']
            self.assertIn('jointReconstruction fluidCubicAdaptive;',(case/'constant/mechanicalProperties').read_text())
            for filename in (('log.first','log.restart') if meta['restart'] else ('log.solver',)):
                log=case/filename;text=log.read_text()
                text+=f'M2A_RECONSTRUCTION scheme=fluidCubicAdaptive moment=1e-15 l1={g["coefficient_l1_max"]} supportMin={g["support_min"]} supportMax={g["support_max"]} expanded={g["expanded_markers"]}\n'
                if filename!='log.restart':text+=f'M2A_SPATIAL_REFERENCE markers={g["markers"]} wallRms={r["wall_rms_m_s"]} wallMax={r["wall_max_m_s"]} stressRatio={r["stress_load_ratio"]}\n'
                log.write_text(text)
            if meta['restart']:
                checkpoint=case/'0.002/couplingState';checkpoint.write_text(checkpoint.read_text()+'\njointReconstruction fluidCubicAdaptive;\n')
            self.assertTrue(summary.audit.evaluate(case)[0]['passed'])
            log=case/('log.restart' if meta['restart'] else 'log.solver');text=log.read_text()
            log.write_text(text.replace('scheme=fluidCubicAdaptive','scheme=trilinear'))
            with self.assertRaisesRegex(ValueError,'reconstruction evidence'):summary.audit.evaluate(case)
            log.write_text(text)
            for bad in (text.replace('moment=1e-15','moment=nan'),re.sub(r'l1=\S+','l1=9',text)):
                log.write_text(bad)
                with self.assertRaisesRegex(ValueError,'reconstruction geometry'):summary.audit.evaluate(case)
            log.write_text(text)
            if meta['restart']:
                checkpoint.write_text(checkpoint.read_text().replace('jointReconstruction fluidCubicAdaptive;','jointReconstruction trilinear;'))
                self.assertFalse(summary.audit.evaluate(case)[0]['checks']['joint_reconstruction_checkpoint'])
            else:
                log.write_text(re.sub(r'wallRms=\S+','wallRms=nan',text))
                with self.assertRaisesRegex(ValueError,'Nonfinite'):summary.audit.evaluate(case)



if __name__=='__main__':unittest.main()
