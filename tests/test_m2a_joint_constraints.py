import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('projection_reference',
    Path(__file__).parents[1]/'scripts/constraint_projection_reference.py')
oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)


class JointConstraintTests(unittest.TestCase):
    def test_joint_solution_and_stationarity(self):
        initial=[1.,0.,0.];mass=[2.,3.,4.]
        rows=[[1.,-1.,0.],[0.,1.,1.]];target=[0.,1.]
        u,mu=oracle.project(initial,mass,rows,target)
        for row,value in zip(rows,target):
            self.assertAlmostEqual(sum(a*b for a,b in zip(row,u)),value,places=13)
        for i in range(3):
            self.assertAlmostEqual(mass[i]*(u[i]-initial[i])+sum(row[i]*x for row,x in zip(rows,mu)),0.)

    def test_sequential_projection_can_break_first_constraint(self):
        initial=[1.,0.,0.];mass=[2.,3.,4.]
        divergence=[1.,-1.,0.];wall=[0.,1.,1.]
        first,_=oracle.project(initial,mass,[divergence],[0.])
        second,_=oracle.project(first,mass,[wall],[1.])
        self.assertGreater(abs(sum(x*y for x,y in zip(divergence,second))),.1)

    def test_adjoint_work_and_force_conservation(self):
        j=[[.25,.75,0.],[0.,.4,.6]];area=[2.,3.];force=[5.,-1.]
        volumes=[.2,.4,.6];u=[1.,2.,-3.]
        spread=oracle.spread(j,area,force,volumes)
        grid_work=sum(a*b*c for a,b,c in zip(u,spread,volumes))
        surface_work=sum(w*f*sum(a*b for a,b in zip(row,u)) for row,w,f in zip(j,area,force))
        self.assertAlmostEqual(grid_work,surface_work)
        self.assertAlmostEqual(sum(f*v for f,v in zip(spread,volumes)),sum(w*f for w,f in zip(area,force)))

    def test_dependent_constraints_rejected(self):
        with self.assertRaisesRegex(ValueError,'Dependent'):
            oracle.project([0.,0.],[1.,1.],[[1.,1.],[2.,2.]],[0.,1.])


if __name__=='__main__':unittest.main()
