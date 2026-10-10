"""Full Cartesian image interpolation on analytic fields; not CFD validation.

Separates interpolation/extension slip from one-sided pressure/shear quadrature.
No pressure projection, solved velocity, DEM evolution or force ledger is present.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import surface_reference as ref


def analytic(r, mode):
    distance=math.sqrt(sum(x*x for x in r)); a=.01/distance
    if mode=='rotate':
        return [-.1*r[1]*a**3,.1*r[0]*a**3,0.,0.]
    n=[x/distance for x in r]
    velocity=[(1-.75*a-.25*a**3)*u+(-.75*a+.75*a**3)*.01*n[0]*v
              for u,v in zip((.01,0.,0.),n)]
    return velocity+[-1.5*.1*.01*.01*r[0]/distance**3]


def sample(point, origin, h, value):
    coordinate=[(x-o)/h for x,o in zip(point,origin)]
    lower=[math.floor(x) for x in coordinate]
    table={tuple(i+j for i,j in zip(lower,offset)):value(tuple(i+j for i,j in zip(lower,offset)))
           for offset in itertools.product((0,1),repeat=3)}
    return ref.trilinear(point,origin,h,table)


def fixture(mode,n,phase,order):
    h=.12/n; radius=.01; origin=[-.06+h/2]*3;centre=[phase*h]*3
    omega=[0.,0.,.1 if mode=='rotate' else 0.]
    fluid_cache={}; ghost_cache={}
    def fluid(index):
        if index not in fluid_cache:
            r=[o+i*h-c for o,i,c in zip(origin,index,centre)]
            if sum(x*x for x in r)<=radius**2:
                raise ValueError('Image stencil contains a solid centre')
            fluid_cache[index]=analytic(r,mode)
        return fluid_cache[index]
    def point_at(normal,distance):
        return [c+distance*x for c,x in zip(centre,normal)]
    def ghost(index):
        if index not in ghost_cache:
            point=[o+i*h for o,i in zip(origin,index)]
            r=[x-c for x,c in zip(point,centre)]
            distance=math.sqrt(sum(x*x for x in r))
            if distance>=radius: result=analytic(r,mode)
            elif distance<=radius-math.sqrt(3)*h or distance<1e-30:
                result=ref.cross(omega,r)+[0.]
            else:
                normal=[x/distance for x in r]
                ub=ref.cross(omega,[radius*x for x in normal])
                a=sample(point_at(normal,radius+2*h),origin,h,fluid)[:3]
                if order=='linear': velocity=ref.ghost_target(ub,a,distance-radius,h)
                else:
                    b=sample(point_at(normal,radius+3*h),origin,h,fluid)[:3]
                    velocity=ref.quadratic_ghost_target(ub,a,b,distance-radius,h)
                result=velocity+[0.]
            ghost_cache[index]=result
        return ghost_cache[index]
    slip=0.;pressure=[0.]*3;viscous=[0.]*3;torque=[0.]*3
    for normal,weight in ref.sphere_quadrature(radius):
        ub=ref.cross(omega,[radius*x for x in normal])
        wall=sample(point_at(normal,radius),origin,h,ghost)[:3]
        slip+=weight*sum((u-v)**2 for u,v in zip(wall,ub))
        a,b,c=[sample(point_at(normal,radius+i*h),origin,h,fluid) for i in (2,3,4)]
        derivative=ref.radial_derivative(ub,a[:3],b[:3],c[:3],h)
        fp,fv=ref.traction(normal,derivative,omega,ref.pressure_wall(a[3],b[3],c[3]),1.,.1)
        pressure=[x+weight*y for x,y in zip(pressure,fp)]
        viscous=[x+weight*y for x,y in zip(viscous,fv)]
        t=ref.cross([radius*x for x in normal],[x+y for x,y in zip(fp,fv)])
        torque=[x+weight*y for x,y in zip(torque,t)]
    force_ref=6*math.pi*.1*radius*.01; torque_ref=-8*math.pi*.1*radius**3*.1
    return {'mode':mode,'mesh_n':n,'phase_h':phase,'order':order,
            'interpolated_image_wall_slip_relative':math.sqrt(slip/(4*math.pi*radius**2))/(.01 if mode=='fixed' else .001),
            'pressure_force':pressure,'viscous_force':viscous,'stress_torque':torque,
            'pressure_drag_ratio':pressure[0]/force_ref,
            'viscous_drag_ratio':viscous[0]/force_ref,
            'stress_load_ratio':(pressure[0]+viscous[0])/force_ref if mode=='fixed' else torque[2]/torque_ref}


def run():
    return {'kind':'analytic_cartesian_wall_operator_audit_only',
            'excluded':['CFD momentum/pressure solution','ghost-pressure coupling','force ledger','native C++ execution'],
            'cases':[fixture(mode,n,phase,order) for mode in ('fixed','rotate')
                     for n in (36,48,64) for phase in (0.,.25) for order in ('linear','quadratic')]}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
