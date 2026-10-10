"""Spatial interpolation on analytic Stokes fields, not a CFD execution.

Images sample the exact continuum field; this omits their Cartesian interpolation
and all fluid-solver errors. Results diagnose geometric/profile bias only.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import surface_reference as reference


def fixture(mode,ngrid,phase):
    h,radius=.12/ngrid,.01
    origin=[-.06+h/2]*3;centre=[phase*h]*3
    def exact(r):
        distance=math.sqrt(sum(x*x for x in r));normal=[x/distance for x in r];a=radius/distance
        if mode=='rotate':return [x*radius*a*a for x in reference.cross([0,0,.1],normal)]
        return [(1-.75*a-.25*a**3)*u+(-.75*a+.75*a**3)*.01*normal[0]*n
                for u,n in zip([.01,0,0],normal)]
    sums=[0.,0.]
    for normal,weight in reference.sphere_quadrature(radius):
        wall=[c+radius*n for c,n in zip(centre,normal)]
        ub=[0.]*3 if mode=='fixed' else [x*radius for x in reference.cross([0,0,.1],normal)]
        lower=[math.floor((x-o)/h) for x,o in zip(wall,origin)]
        tables=[{},{}]
        for offsets in itertools.product((0,1),repeat=3):
            index=tuple(i+j for i,j in zip(lower,offsets))
            point=[o+i*h for o,i in zip(origin,index)]
            r=[x-c for x,c in zip(point,centre)];distance=math.sqrt(sum(x*x for x in r))
            if distance>=radius:
                values=[exact(r)]*2
            elif distance<=radius-math.sqrt(3)*h or distance==0:
                rigid=[0.]*3 if mode=='fixed' else reference.cross([0,0,.1],r)
                values=[rigid]*2
            else:
                n=[x/distance for x in r]
                b=[0.]*3 if mode=='fixed' else [x*radius for x in reference.cross([0,0,.1],n)]
                a=exact([(radius+2*h)*x for x in n]);c=exact([(radius+3*h)*x for x in n])
                values=[reference.ghost_target(b,a,distance-radius,h),
                        reference.quadratic_ghost_target(b,a,c,distance-radius,h)]
            for order in range(2):tables[order][index]=values[order]
        for order in range(2):
            value=reference.trilinear(wall,origin,h,tables[order])
            sums[order]+=weight*sum((a-b)**2 for a,b in zip(value,ub))
    scale=.01 if mode=='fixed' else .001
    return {'mode':mode,'mesh_n':ngrid,'phase_h':phase,
            'linear_relative_slip':math.sqrt(sums[0]/(4*math.pi*radius**2))/scale,
            'quadratic_relative_slip':math.sqrt(sums[1]/(4*math.pi*radius**2))/scale}


def run():
    return {'kind':'analytic_field_geometry_diagnostic_only',
            'not_included':['CFD errors','Cartesian image interpolation errors','force validation'],
            'cases':[fixture(mode,n,phase) for mode in ('fixed','rotate')
                     for n in (36,48,64) for phase in (0.,.25)]}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for case in result['cases']:print(case)
