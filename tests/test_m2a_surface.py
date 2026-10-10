"""Surface formulas and independent gates. Synthetic ledgers are not CFD proof."""
import csv
import importlib.util
import itertools
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
import surface_reference as reference


def load(name,filename):
    spec=importlib.util.spec_from_file_location(name,SCRIPTS/filename)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


prepare=load('prepare_surface_test','prepare-m2a-surface.py')
audit=load('audit_surface_test','summarize-m2a-surface.py')


class FormulaTests(unittest.TestCase):
    def test_affine_interpolation_and_partition(self):
        values={i:[1+2*i[0]-3*i[1]+4*i[2],1.] for i in itertools.product((0,1),repeat=3)}
        self.assertEqual(reference.trilinear([.2,.4,.6],[0]*3,1,values),[2.6,1.])
        self.assertEqual(reference.trilinear([0,0,0],[0]*3,1,values),values[(0,0,0)])

    def test_image_stencil_is_fluid_side(self):
        radius=.01
        for ngrid in (36,48,64):
            h=.12/ngrid;origin=-.06+h/2
            for phase in (0.,.25):
                centre=[phase*h]*3
                for normal,weight in reference.sphere_quadrature(radius):
                    image=[c+(radius+2*h)*n for c,n in zip(centre,normal)]
                    lower=[math.floor((x-origin)/h) for x in image]
                    for offsets in itertools.product((0,1),repeat=3):
                        corner=[origin+(i+j)*h for i,j in zip(lower,offsets)]
                        self.assertGreater(sum((x-c)**2 for x,c in zip(corner,centre)),radius**2)

    def test_ghost_linear_profile_and_radial_cubic(self):
        h=.003
        wall=[.2,-.1,.3]; slope=[2.,3.,-4.]
        image=[a+2*h*b for a,b in zip(wall,slope)]
        target=reference.ghost_target(wall,image,-.75*h,h)
        for a,b,c in zip(target,wall,slope):
            self.assertAlmostEqual(a,b-.75*h*c)
        def profile(t):
            return [a+b*t+7*t*t-13*t**3 for a,b in zip(wall,slope)]
        derivative=reference.radial_derivative(wall,profile(2*h),profile(3*h),profile(4*h),h)
        for a,b in zip(derivative,slope):self.assertAlmostEqual(a,b,places=11)
        self.assertAlmostEqual(reference.pressure_wall(1+2*(2*h)+3*(2*h)**2,
            1+2*(3*h)+3*(3*h)**2,1+2*(4*h)+3*(4*h)**2),1)

    def test_quadratic_ghost_profile_and_shared_wall(self):
        h=.003
        def profile(t):return [.2+2*t+5*t*t,-.1-4*t-7*t*t,.3+t]
        for distance in (-math.sqrt(3)*h,-.5*h,0.):
            actual=reference.quadratic_ghost_target(profile(0),profile(2*h),profile(3*h),distance,h)
            for a,b in zip(actual,profile(distance)):self.assertAlmostEqual(a,b,places=13)

    def test_quadrature_closure_and_pressure_gauge(self):
        nodes=list(reference.sphere_quadrature(.01))
        self.assertEqual(len(nodes),288)
        area=4*math.pi*.01**2
        self.assertAlmostEqual(sum(w for n,w in nodes)/area,1,places=13)
        for axis in range(3):
            self.assertLess(abs(sum(n[axis]*w for n,w in nodes)/area),1e-14)
            self.assertAlmostEqual(sum(n[axis]**2*w for n,w in nodes)/area,1/3,places=13)
        force=[sum(-103*rho*n[a]*w for n,w in nodes) for a in range(3) for rho in [2.]]
        self.assertLess(math.sqrt(sum(x*x for x in force)),1e-15)

    def test_stokes_pressure_viscosity_and_rotating_torque(self):
        r,rho,nu=.01,2.,.1
        u=[.01,0,0];omega=[0,0,.1]
        fp,fv,torque=[0.]*3,[0.]*3,[0.]*3
        for n,w in reference.sphere_quadrature(r):
            dot=sum(a*b for a,b in zip(u,n))
            radial=[1.5/r*(a-dot*b) for a,b in zip(u,n)]
            p,v=reference.traction(n,radial,[0]*3,-1.5*nu*dot/r,rho,nu)
            fp=[a+b*w for a,b in zip(fp,p)];fv=[a+b*w for a,b in zip(fv,v)]
            radial=[-2*x for x in reference.cross(omega,n)]
            p,v=reference.traction(n,radial,omega,0,rho,nu)
            t=reference.cross([r*x for x in n],v)
            torque=[a+b*w for a,b in zip(torque,t)]
        self.assertAlmostEqual(fp[0]/(2*math.pi*rho*nu*r*u[0]),1,places=13)
        self.assertAlmostEqual(fv[0]/(4*math.pi*rho*nu*r*u[0]),1,places=13)
        self.assertAlmostEqual(-torque[2]/(8*math.pi*rho*nu*r**3*omega[2]),1,places=13)
        _,rigid=reference.traction([1,0,0],reference.cross(omega,[1,0,0]),omega,0,rho,nu)
        self.assertEqual(rigid,[0,0,0])


class SurfaceGateTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()

    def fixture(self,name,ratio=1.):
        case=self.root/name;meta=prepare.prepare(case,name)
        (case/'result.txt').write_text('exit_status=0\n')
        for log in ('log.first','log.restart') if meta['restart'] else ('log.solver',):
            (case/log).write_text('M2A_CONSTRAINT scheme=surfaceExtension\nM2A_RECONSTRUCTION order=quadratic\nM2A_EXECUTION_COMPLETE\n')
        if meta['restart']:
            (case/'0.002').mkdir()
            for filename in ('couplingState','dem.restart'):(case/'0.002'/filename).write_text('synthetic')
        r=meta['radius']
        exact=(6*math.pi*meta['rho']*meta['nu']*r*.01 if meta['mode']=='fixed'
               else 8*math.pi*meta['rho']*meta['nu']*r**3*.1)
        force=[exact*ratio,0,0] if meta['mode']=='fixed' else [0.,0,0]
        torque=[0.,0,-exact*ratio] if meta['mode']=='rotate' else [0.,0,0]
        rows=[]
        for window in range(1,41):
            row={'time':window*1e-4,'id':7,'owner':0,'volume_error':.005,'slip_rms':0.,
                 'div_max':0.,'clock_error':0.,'covered_ranks':meta['nprocs'],
                 'pressure_correctors':8,'momentum_equation_impulse_L1':0.,'support_volume_ratio':1.,
                 'surface_correctors':2,'surface_target_defect':1e-8,'pressure_correctors_total':16,
                 'window_wall_seconds':.1}
            for a,c in enumerate('xyz'):
                row[c]=meta['initial_centre'][a];row['v'+c]=meta['initial_velocity'][a]
                row['w'+c]=meta['initial_omega'][a];row['f'+c]=force[a];row['t'+c]=torque[a]
                for prefix in ('fluid_p','old_fluid_p','inertia_f','inertia_t'):row[prefix+c]=0.
                row['boundary_f'+c]=force[a];row['constraint_f'+c]=-force[a];row['constraint_t'+c]=-torque[a]
                row['stress_pressure_f'+c]=force[a]/3;row['stress_viscous_f'+c]=force[a]*2/3
                row['stress_pressure_t'+c]=0.;row['stress_viscous_t'+c]=torque[a]
            rows.append(row)
        self.write(case,rows)
        return case,rows,meta

    def write(self,case,rows):
        with (case/'mechanical-history.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)

    def test_matrix_and_phase_positions(self):
        for name in prepare.CASES:
            case,rows,meta=self.fixture(name)
            self.assertIn(f'({meta["mesh_n"]} {meta["mesh_n"]} {meta["mesh_n"]})',
                          (case/'system/blockMeshDict').read_text())
            centre=' '.join(format(x,'.17g') for x in meta['initial_centre'])
            self.assertIn('7 1 0.02 1000 '+centre,(case/'particle.data').read_text())
        self.assertTrue(audit.summarize(self.root)['passed'])

    def test_wall_slip_target_and_budget_failures(self):
        case,rows,meta=self.fixture('fixed-coarse')
        rows[0]['slip_rms']=.001;rows[0]['surface_target_defect']=2e-7
        rows[0]['surface_correctors']=33;self.write(case,rows)
        result,_,_=audit.evaluate(case)
        for key in ('slip','surface_target_convergence','surface_iteration_budget'):
            self.assertFalse(result['checks'][key])

    def test_native_scheme_marker_required(self):
        case,rows,meta=self.fixture('fixed-coarse')
        (case/'log.solver').write_text('M2A_EXECUTION_COMPLETE\n')
        with self.assertRaisesRegex(ValueError,'scheme marker'):
            audit.evaluate(case)

    def test_continuity_and_reconstruction_cannot_be_hidden(self):
        case,rows,meta=self.fixture('fixed-coarse')
        rows[0]['div_max']=2e-7;self.write(case,rows)
        result,_,_=audit.evaluate(case)
        self.assertTrue(result['checks']['divergence'])
        self.assertFalse(result['checks']['continuity_convergence'])
        (case/'log.solver').write_text('M2A_CONSTRAINT scheme=surfaceExtension\nM2A_EXECUTION_COMPLETE\n')
        with self.assertRaisesRegex(ValueError,'quadratic'):audit.evaluate(case)

    def test_fine_force_stress_and_refinement_failure(self):
        for name in prepare.CASES:self.fixture(name,1.12 if 'finer' in name else 1.)
        case=self.root/'fixed-finer'
        entry,rows,meta=audit.evaluate(case)
        self.assertFalse(entry['checks']['fine64_reference'])
        rows[-1]['stress_viscous_fx']+=20*.12*6*math.pi*.1*.01*.01
        self.write(case,rows)
        self.assertFalse(audit.evaluate(case)[0]['checks']['fine64_stress_ledger'])
        summary=audit.summarize(self.root)
        self.assertFalse(summary['comparisons']['fixed-fine-to-finer']['passed'])

    def test_phase_and_missing_comparison_remain_failed(self):
        for name in prepare.CASES:
            if name not in ('fixed-restart','fixed-mpi2'):
                self.fixture(name,1.04 if name.endswith('-offset') else 1.)
        summary=audit.summarize(self.root)
        self.assertFalse(summary['comparisons']['fixed-finer-phase']['passed'])
        self.assertEqual(summary['comparisons']['fixed-mpi2']['status'],'not_evaluated')
        self.assertFalse(summary['passed'])


if __name__=='__main__':unittest.main()
