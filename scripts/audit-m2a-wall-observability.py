"""Geometric null-direction counterexample for fluid-only J and wall sampling Q.

No FV pressure projection or momentum solve: this is a compatibility bound,
not proof that the exhibited velocity occurs in the solved CFD trajectory.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import surface_reference as ref


def fixture(n,phase):
    h=.12/n;radius=.01;origin=[-.06+h/2]*3;centre=[phase*h]*3
    total=interior=0.;example=None;touches=0
    for normal,area in ref.sphere_quadrature(radius):
        point=[c+radius*x for c,x in zip(centre,normal)]
        coordinates=[(x-o)/h for x,o in zip(point,origin)]
        lower=[math.floor(x) for x in coordinates]
        fractions=[x-i for x,i in zip(coordinates,lower)]
        energy=inside=0.
        for delta in itertools.product((0,1),repeat=3):
            index=tuple(i+d for i,d in zip(lower,delta))
            weight=math.prod(f if d else 1-f for f,d in zip(fractions,delta))
            energy+=weight*weight
            position=[o+i*h-c for o,i,c in zip(origin,index,centre)]
            if sum(x*x for x in position)<=radius**2:
                inside+=weight*weight
                if example is None or weight>example['independent_wall_response']:
                    example={'interior_cell':list(index),'cell_radius_m':math.sqrt(sum(x*x for x in position)),
                             'unit_cell_velocity_m_s':1.,'fluid_only_J_response_m_s':0.,
                             'independent_wall_response':weight}
        if inside>0:touches+=1
        total+=area*energy;interior+=area*inside
    return {'mesh_n':n,'phase_h':phase,'quadrature_nodes':288,'nodes_touching_interior':touches,
            'interior_fraction_of_area_weighted_Q_frobenius_energy':interior/total,
            'lower_bound_Q_nullJ_frobenius_over_Q':math.sqrt(interior/total),
            'null_direction_example':example,
            'interpretation':'Interior columns of fluid-only J are zero, but Q has nonzero interior columns. J U=Ub alone cannot guarantee independent Q wall values. Bound uses unweighted cell velocity norm.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result={'scope':__doc__,'cases':[fixture(n,p) for n in (36,48,64) for p in (0.,.25)]}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
