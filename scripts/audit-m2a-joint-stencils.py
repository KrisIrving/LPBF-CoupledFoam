"""Independent geometric J/S/rank audit; no native CFD execution."""
import argparse
import itertools
import json
import math
from pathlib import Path


def fixture(n,phase):
    h=.12/n;radius=.01;centre=[phase*h]*3;origin=[-.06+h/2]*3
    count=math.ceil(4*math.pi*radius**2/(2.25*h*h))
    rows=[];points=[];partition=0.;coordinate=0.
    for i in range(count):
        z=1-2*(i+.5)/count;phi=math.pi*(3-math.sqrt(5))*i
        point=[centre[0]+radius*math.sqrt(1-z*z)*math.cos(phi),
               centre[1]+radius*math.sqrt(1-z*z)*math.sin(phi),centre[2]+radius*z]
        x=[(p-o)/h for p,o in zip(point,origin)]
        lower=[math.floor(a) for a in x];fraction=[a-b for a,b in zip(x,lower)]
        row={tuple(i+j for i,j in zip(lower,offset)):
             math.prod(f if j else 1-f for f,j in zip(fraction,offset))
             for offset in itertools.product((0,1),repeat=3)}
        if any(any(i<0 or i>=n for i in index) for index in row):raise ValueError('Grid bounds')
        partition=max(partition,abs(sum(row.values())-1))
        coordinate=max(coordinate,max(abs(sum(w*(origin[a]+index[a]*h) for index,w in row.items())-point[a]) for a in range(3)))
        points.append(point);rows.append(row)
    diagonal=[sum(w*w for w in row.values()) for row in rows]
    cholesky=[[0.]*count for _ in rows];min_pivot=1.
    for i,row in enumerate(rows):
        for j in range(i+1):
            value=sum(w*rows[j].get(index,0) for index,w in row.items())/math.sqrt(diagonal[i]*diagonal[j])
            value-=sum(cholesky[i][k]*cholesky[j][k] for k in range(j))
            if i==j:
                if value<=1e-8:raise ValueError('Dependent interpolation rows')
                min_pivot=min(min_pivot,value);cholesky[i][j]=math.sqrt(value)
            else:cholesky[i][j]=value/cholesky[j][j]
    # Arbitrary integrated forces test both zero/first moment and work exchange.
    forces=[math.sin(i+.7)*1e-4 for i in range(count)]
    grid={}
    for force,row in zip(forces,rows):
        for index,w in row.items():grid[index]=grid.get(index,0)+w*force
    marker_force=sum(forces);grid_force=sum(grid.values())
    first_moment=max(abs(sum(force*point[a] for force,point in zip(forces,points))-
                         sum(force*(origin[a]+index[a]*h) for index,force in grid.items())) for a in range(3))
    def velocity(index):return math.cos(index[0]+2*index[1]-.5*index[2])
    grid_work=sum(force*velocity(index) for index,force in grid.items())
    marker_work=sum(force*sum(w*velocity(index) for index,w in row.items()) for force,row in zip(forces,rows))
    return {'mesh_n':n,'phase_h':phase,'markers':count,'min_normalised_cholesky_pivot':min_pivot,
            'partition_error':partition,'coordinate_error_m':coordinate,
            'force_exchange_error_N':abs(marker_force-grid_force),
            'first_moment_exchange_error_N_m':first_moment,'work_exchange_error_W':abs(marker_work-grid_work)}


def run():
    cases=[fixture(n,phase) for n in (36,48,64) for phase in (0.,.25)]
    passed=all(c['partition_error']<1e-12 and c['coordinate_error_m']<1e-14
               and c['force_exchange_error_N']<1e-12 and c['first_moment_exchange_error_N_m']<1e-14
               and c['work_exchange_error_W']<1e-14 for c in cases)
    return {'kind':'geometric_stencil_oracle_only','pressure_projected_schur_rank':'not_evaluated',
            'passed':passed,'cases':cases}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for case in result['cases']:print(case)
    raise SystemExit(0 if result['passed'] else 1)
