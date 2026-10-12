"""Independent SVD oracle for complete wall space and modal load transforms.

NumPy is needed only to regenerate this development artifact; case generation
reads the saved six-fixture data. No CFD pressure response is represented.
"""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import surface_reference as ref

spec=importlib.util.spec_from_file_location('complete_sampling',Path(__file__).with_name('audit-m2a-complete-wall.py'))
sampling=importlib.util.module_from_spec(spec);spec.loader.exec_module(sampling)
spec=importlib.util.spec_from_file_location('analytic',Path(__file__).with_name('audit-m2a-wall-operators.py'))
analytic=importlib.util.module_from_spec(spec);spec.loader.exec_module(analytic)


def pivoted_mgs(matrix):
    residual=matrix.copy();basis=[]
    for _ in range(len(matrix)):
        pivot=int(np.argmax(np.sum(residual*residual,axis=1)));value=residual[pivot].copy()
        for _ in range(2):
            for row in basis:value-=np.dot(value,row)*row
        length=float(np.linalg.norm(value))
        if length<=1e-11:break
        value/=length;basis.append(value)
        for _ in range(2):residual-=np.outer(residual@value,value)
    return np.array(basis)


def fixture(n,phase):
    h=.12/n;radius=.01;origin=[-.06+h/2]*3;centre=[phase*h]*3
    points=np.array([[c+radius*x for c,x in zip(centre,normal)] for normal,_ in ref.sphere_quadrature(radius)])
    rows=[sampling.sampling_row(point,origin,h) for point in points]
    cells=sorted(set().union(*(row.keys() for row in rows)))
    q=np.array([[row.get(k,0.) for k in cells] for row in rows])
    positions=np.array([[o+i*h for o,i in zip(origin,index)] for index in cells])
    _,singular,right=np.linalg.svd(q,full_matrices=False)
    rank=int(np.count_nonzero(singular>1e-11));basis=right[:rank]
    mgs=pivoted_mgs(q)
    mgs_coverage=float(np.max(np.linalg.norm(q-(q@mgs.T)@mgs,axis=1)/np.linalg.norm(q,axis=1)))
    mgs_orthogonality=float(np.max(np.abs(mgs@mgs.T-np.eye(len(mgs)))))
    coordinates=q@basis.T
    coverage=float(np.max(np.linalg.norm(q-coordinates@basis,axis=1)/np.linalg.norm(q,axis=1)))
    orthogonality=float(np.max(np.abs(basis@basis.T-np.eye(rank))))
    force=np.array([[1e-6*math.sin(.73*(3*i+a+1)) for a in range(3)] for i in range(rank)])
    cell_force=basis.T@force;sums=basis.sum(axis=1);first=basis@(positions-np.array(centre))
    modal_force=(sums[:,None]*force).sum(axis=0)
    modal_torque=np.cross(first,force).sum(axis=0)
    velocity=np.array([[math.cos(.31*(i+a+1)) for a in range(3)] for i in range(len(cells))])
    errors={'force_N':float(np.max(np.abs(modal_force-cell_force.sum(axis=0)))),
            'torque_N_m':float(np.max(np.abs(modal_torque-np.cross(positions-np.array(centre),cell_force).sum(axis=0)))),
            'work_W':float(abs(np.sum(force*(basis@velocity))-np.sum(cell_force*velocity)))}
    references={}
    for mode in ('fixed','rotate'):
        omega=np.array([0.,0.,.1 if mode=='rotate' else 0.])
        target=np.cross(np.broadcast_to(omega,positions.shape),positions-np.array(centre))
        modal_target=np.cross(np.broadcast_to(omega,first.shape),first)
        rigid_error=float(np.max(np.abs(basis@target-modal_target)))
        values=[]
        for position,r in zip(positions,positions-np.array(centre)):
            values.append(np.cross(omega,r) if np.linalg.norm(r)<=radius else analytic.analytic(r,mode)[:3])
        wall=q@np.array(values);desired=np.cross(np.broadcast_to(omega,points.shape),points-np.array(centre))
        norms=np.linalg.norm(wall-desired,axis=1)
        references[mode]={'wall_rms_m_s':float(np.sqrt(np.mean(norms**2))),
                          'wall_max_m_s':float(np.max(norms)), 'rigid_modal_error_m_s':rigid_error}
    return {'mesh_n':n,'phase_h':phase,'full_wall_nodes':288,'support_cells':len(cells),'modes':rank,
        'smallest_retained_singular':float(singular[rank-1]),
        'largest_discarded_singular':float(singular[rank]) if rank<len(singular) else 0.,
        'coverage_relative':coverage,'orthogonality':orthogonality,'exchange_errors':errors,'references':references,
        'pivoted_mgs_offline':{'modes':len(mgs),'coverage_relative':mgs_coverage,'orthogonality':mgs_orthogonality},
        'admitted_geometry':coverage<=1e-10 and orthogonality<=1e-10
            and len(mgs)==rank and mgs_coverage<=1e-10 and mgs_orthogonality<=1e-10
            and errors['force_N']<=1e-12 and errors['torque_N_m']<=1e-14 and errors['work_W']<=1e-14
            and all(r['rigid_modal_error_m_s']<=1e-10 for r in references.values())}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();cases=[fixture(n,p) for n in (36,48,64) for p in (0.,.25)]
    args.output.write_text(json.dumps({'scope':__doc__,'cases':cases,'all_geometry_admitted':all(c['admitted_geometry'] for c in cases)},indent=2,allow_nan=False)+'\n')
