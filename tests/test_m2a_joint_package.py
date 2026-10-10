"""Generator/gate regressions using synthetic histories, not CFD evidence."""
import importlib.util
import json
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
synthetic=load('joint_legacy_fixture',root/'tests/test_m2a_surface.py')


class JointPackageTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.writer=synthetic.SurfaceGateTests();self.writer.root=self.root
    def tearDown(self):self.temp.cleanup()
    def fixture(self,name):
        case,rows,_=self.writer.fixture(name)
        meta=prepare.prepare(case,name)
        for log in ('log.first','log.restart') if meta['restart'] else ('log.solver',):
            (case/log).write_text('M2A_CONSTRAINT scheme=jointSurface\nM2A_JOINT markers='+str(meta['joint_controls']['markers'])+' minPivot=.9\n'
                'M2A_BOUNDARY treatment=compatibleGauss rawFlux=-4e-14 gaussFlux=7e-18 finalFlux=1e-22 normalCorrection=1e-16\nM2A_EXECUTION_COMPLETE\n')
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
            self.assertNotIn('constraintScheme surfaceExtension;',properties)
            self.assertEqual(meta['package'],'M2A-02C')
        self.assertTrue(summary.summarize(self.root)['passed'])

    def test_legacy_runs_cannot_masquerade_as_joint(self):
        for name in prepare.CASES:self.writer.fixture(name)
        result=summary.summarize(self.root)
        self.assertFalse(result['passed'])
        self.assertTrue(all(not entry['passed'] for entry in result['cases'].values()))

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


if __name__=='__main__':unittest.main()
