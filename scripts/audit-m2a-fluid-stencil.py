"""Fluid-side polynomial J/transpose-S preflight; not projected CFD evidence."""
import argparse
import itertools
import json
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import surface_reference as ref
from importlib.util import spec_from_file_location,module_from_spec
spec=spec_from_file_location('analytic_wall',Path(__file__).with_name('audit-m2a-wall-operators.py'))
wall=module_from_spec(spec);spec.loader.exec_module(wall)
spec=spec_from_file_location('baseline_spatial',Path(__file__).with_name('audit-m2a-spatial-consistency.py'))
baseline=module_from_spec(spec);spec.loader.exec_module(baseline)


def weights(offsets,degree=2):
    """Weighted polynomial evaluation at origin, twice-MGS QR; no regularization."""
    if degree not in (2,3):raise ValueError('Degree must be 2 or 3')
    dimension=10 if degree==2 else 20
    basis=[[1.,x,y,z,x*x,x*y,x*z,y*y,y*z,z*z,x*x*x,x*x*y,x*x*z,x*y*y,x*y*z,x*z*z,y*y*y,y*y*z,y*z*z,z*z*z][:dimension] for x,y,z in offsets]
    d=[math.exp(-sum(x*x for x in r)/2) for r in offsets]
    q=[];r=[[0.]*dimension for _ in range(dimension)]
    for j in range(dimension):
        v=[math.sqrt(w)*b[j] for w,b in zip(d,basis)]
        original=math.sqrt(sum(x*x for x in v))
        for _ in range(2):
            for i,column in enumerate(q):
                projection=sum(a*b for a,b in zip(column,v));r[i][j]+=projection
                v=[a-projection*b for a,b in zip(v,column)]
        r[j][j]=math.sqrt(sum(x*x for x in v))
        if r[j][j]<=max(original,1.)*1e-10:raise ValueError(f'Deficient fluid polynomial stencil at basis column {j}; pivot={r[j][j]}')
        q.append([x/r[j][j] for x in v])
    z=[]
    for i in range(dimension):z.append(((1. if i==0 else 0.)-sum(r[j][i]*z[j] for j in range(i)))/r[i][i])
    result=[math.sqrt(d[i])*sum(q[j][i]*z[j] for j in range(dimension)) for i in range(len(offsets))]
    defect=max(abs(sum(w*b[j] for w,b in zip(result,basis))-(1. if j==0 else 0.)) for j in range(dimension))
    if defect>1e-10:raise ValueError('Polynomial reproduction failed')
    return result,defect


def fixture(n,phase,degree=2,adaptive=False):
    h=.12/n;radius=.01;centre=[phase*h]*3;origin=[-.06+h/2]*3
    count=math.ceil(4*math.pi*radius**2/(2.25*h*h));rows=[];points=[];moment=0.;l1=0.;fixed=[];rotate=[];expanded=0
    for i in range(count):
        z=1-2*(i+.5)/count;phi=math.pi*(3-math.sqrt(5))*i
        point=[c+radius*x for c,x in zip(centre,[math.sqrt(1-z*z)*math.cos(phi),math.sqrt(1-z*z)*math.sin(phi),z])]
        lower=[math.floor((x-o)/h) for x,o in zip(point,origin)]
        for support in ((2.5,3.5) if adaptive else (2.5,)):
            indices=[];offsets=[];width=math.ceil(support)
            for delta in itertools.product(range(1-width,width+1),repeat=3):
                index=tuple(i+j for i,j in zip(lower,delta));cell=[o+i*h for o,i in zip(origin,index)]
                offset=[(x-p)/h for x,p in zip(cell,point)]
                if sum(x*x for x in offset)>support**2 or sum((x-c)**2 for x,c in zip(cell,centre))<=radius**2:continue
                if any(i<0 or i>=n for i in index):raise ValueError('Grid bounds')
                indices.append(index);offsets.append(offset)
            try:
                w,error=weights(offsets,degree)
                if sum(abs(v) for v in w)>8:raise ValueError('Unbounded reconstruction')
                break
            except ValueError as error:
                if not adaptive or support==3.5:
                    levels=[len({index[a] for index in indices}) for a in range(3)]
                    raise ValueError(f'marker={i}, support={len(indices)}, axis_levels={levels}: {error}') from error
        moment=max(moment,error);l1=max(l1,sum(abs(x) for x in w))
        if support>2.5:expanded+=1
        row=dict(zip(indices,w));rows.append(row);points.append(point)
        for mode,errors in (('fixed',fixed),('rotate',rotate)):
            actual=[sum(v*wall.analytic([o+i*h-c for o,i,c in zip(origin,index,centre)],mode)[a] for index,v in row.items()) for a in range(3)]
            desired=[0.,0.,0.] if mode=='fixed' else ref.cross([0,0,.1],[p-c for p,c in zip(point,centre)])
            errors.append(math.sqrt(sum((a-b)**2 for a,b in zip(actual,desired))))
    diagonal=[sum(w*w for w in row.values()) for row in rows];factor=[[0.]*count for _ in rows];pivot=1.
    for i,row in enumerate(rows):
        for j in range(i+1):
            value=sum(w*rows[j].get(index,0) for index,w in row.items())/math.sqrt(diagonal[i]*diagonal[j])
            value-=sum(factor[i][k]*factor[j][k] for k in range(j))
            if i==j:
                if value<=1e-8:raise ValueError('Dependent fluid-side J rows')
                pivot=min(pivot,value);factor[i][j]=math.sqrt(value)
            else:factor[i][j]=value/factor[j][j]
    forces=[math.sin(i+.7)*1e-4 for i in range(count)];grid={}
    for force,row in zip(forces,rows):
        for index,w in row.items():grid[index]=grid.get(index,0)+w*force
    force_error=abs(sum(forces)-sum(grid.values()))
    first=max(abs(sum(f*p[a] for f,p in zip(forces,points))-sum(f*(origin[a]+index[a]*h) for index,f in grid.items())) for a in range(3))
    field=lambda index:math.cos(index[0]+2*index[1]-.5*index[2])
    work=abs(sum(f*field(index) for index,f in grid.items())-sum(f*sum(w*field(index) for index,w in row.items()) for f,row in zip(forces,rows)))
    return {'mesh_n':n,'phase_h':phase,'markers':count,'expanded_markers':expanded,'support_min':min(map(len,rows)),
        'support_max':max(map(len,rows)),'coefficient_l1_max':l1,'polynomial_moment_defect':moment,
        'min_normalised_cholesky_pivot':pivot,'fixed_wall_rms_m_s':math.sqrt(sum(x*x for x in fixed)/count),
        'rotate_wall_rms_m_s':math.sqrt(sum(x*x for x in rotate)/count),
        'fixed_wall_max_m_s':max(fixed),'rotate_wall_max_m_s':max(rotate),
        'force_exchange_error_N':force_error,'first_moment_exchange_error_N_m':first,'work_exchange_error_W':work}


def run():
    candidates=[]
    for degree,adaptive in ((2,False),(3,False),(3,True)):
        cases=[]
        for n in (36,48,64):
            for phase in (0.,.25):
                try:
                    case=fixture(n,phase,degree,adaptive)
                    checks={'moments':case['polynomial_moment_defect']<=1e-10,
                        'bounded_coefficients':case['coefficient_l1_max']<=8,
                        'force_exchange':case['force_exchange_error_N']<=1e-12,
                        'first_moment_exchange':case['first_moment_exchange_error_N_m']<=1e-14,
                        'work_exchange':case['work_exchange_error_W']<=1e-14}
                    for mode in ('fixed','rotate'):
                        old=baseline.fixture(mode,n,phase)['wall_rms_m_s']
                        case[mode+'_baseline_wall_rms_m_s']=old
                        checks[mode+'_wall_improvement']=case[mode+'_wall_rms_m_s']<old
                    case.update(checks=checks,admitted=all(checks.values()))
                except ValueError as error:case={'mesh_n':n,'phase_h':phase,'admitted':False,'error':str(error)}
                cases.append(case)
        candidates.append({'degree':degree,'adaptive':adaptive,'cases':cases,'admitted_to_native':all(c['admitted'] for c in cases)})
    return {'scope':__doc__,'pressure_projected_schur_rank':'not_evaluated',
            'support_radius_h':2.5,'weight':'exp(-distance_h^2/2)',
            'admission':'All six geometries improve both analytic wall RMS; moment/exchange/coefficient gates. Not physical validation.',
            'candidates':candidates}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=run()
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')

