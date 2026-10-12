"""Preflight full independent wall quadrature as a paired constraint space.

Pivoted Cholesky selects original sampling rows, without regularization.
This checks unprojected geometry only, not the FV pressure Schur operator.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import surface_reference as ref


def sampling_row(point,origin,h):
    coordinates=[(p-o)/h for p,o in zip(point,origin)]
    lower=[math.floor(x) for x in coordinates]
    fractions=[x-i for x,i in zip(coordinates,lower)]
    return {tuple(i+d for i,d in zip(lower,delta)):
            math.prod(f if d else 1-f for f,d in zip(fractions,delta))
            for delta in itertools.product((0,1),repeat=3)}


def select_rows(rows,pivot_limit=1e-8):
    norms=[math.sqrt(sum(v*v for v in row.values())) for row in rows]
    gram=[[sum(v*other.get(k,0.) for k,v in row.items())/(norms[i]*norms[j])
           for j,other in enumerate(rows)] for i,row in enumerate(rows)]
    diagonal=[gram[i][i] for i in range(len(rows))]
    factors=[[] for _ in rows];chosen=[];pivots=[]
    while len(chosen)<len(rows):
        candidates=[i for i in range(len(rows)) if i not in chosen]
        index=max(candidates,key=lambda i:diagonal[i]);pivot=diagonal[index]
        if pivot<=pivot_limit:break
        chosen.append(index);pivots.append(pivot)
        column=[(gram[i][index]-sum(a*b for a,b in zip(factors[i],factors[index])))/math.sqrt(pivot)
                for i in range(len(rows))]
        for i,value in enumerate(column):
            factors[i].append(value);diagonal[i]-=value*value
    # Do not use sqrt(Gram diagonal subtraction) as the final span error:
    # cancellation at machine precision would manufacture ~1e-8 errors.
    maximum=0.
    for row_index,row in enumerate(rows):
        y=[]
        for i,index in enumerate(chosen):
            y.append((gram[row_index][index]-sum(factors[index][k]*y[k] for k in range(i)))/factors[index][i])
        coefficients=[0.]*len(chosen)
        for i in range(len(chosen)-1,-1,-1):
            coefficients[i]=(y[i]-sum(factors[chosen[j]][i]*coefficients[j]
                                     for j in range(i+1,len(chosen))))/factors[chosen[i]][i]
        residual={k:v/norms[row_index] for k,v in row.items()}
        for coefficient,index in zip(coefficients,chosen):
            for k,v in rows[index].items():residual[k]=residual.get(k,0.)-coefficient*v/norms[index]
        maximum=max(maximum,math.sqrt(sum(v*v for v in residual.values())))
    return chosen,pivots,maximum


def fixture(n,phase):
    h=.12/n;radius=.01;origin=[-.06+h/2]*3;centre=[phase*h]*3
    points=[[c+radius*x for c,x in zip(centre,normal)] for normal,_ in ref.sphere_quadrature(radius)]
    rows=[sampling_row(point,origin,h) for point in points]
    chosen,pivots,remainder=select_rows(rows)
    moment=max(abs(sum(row.values())-1) for row in rows)
    first=max(abs(sum(w*(origin[a]+index[a]*h) for index,w in row.items())-point[a])/h
              for row,point in zip(rows,points) for a in range(3))
    return {'mesh_n':n,'phase_h':phase,'full_wall_nodes':len(rows),'independent_rows':len(chosen),
        'selected_indices':chosen,'minimum_selected_pivot':min(pivots),
        'max_unresolved_row_relative_norm':remainder,'constant_defect':moment,'first_moment_defect_h':first,
        'covers_full_Q_to_1e_minus_10':remainder<=1e-10,
        'admitted_to_native':remainder<=1e-10 and moment<=1e-10 and first<=1e-10,
        'scope':'Selected J=original Q rows, S=transpose/rhoV. Independence is unprojected only; actual FV/pressure/traction not tested.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result={'scope':__doc__,'pivot_limit':1e-8,
        'cases':[fixture(n,p) for n in (36,48,64) for p in (0.,.25)]}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
