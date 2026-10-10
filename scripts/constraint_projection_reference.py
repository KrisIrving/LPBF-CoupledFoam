"""Small dense algebra oracle for joint constraints; NOT an OpenFOAM solver.

min 1/2 (u-u*)^T M (u-u*) subject to C u = b.
Production must use the actual FV face-flux pressure operator, not this toy C.
"""


def solve(matrix, rhs):
    n=len(rhs)
    augmented=[list(row)+[value] for row,value in zip(matrix,rhs)]
    for col in range(n):
        pivot=max(range(col,n),key=lambda row:abs(augmented[row][col]))
        if abs(augmented[pivot][col])<1e-12:
            raise ValueError('Dependent or ill-conditioned constraint rows')
        augmented[col],augmented[pivot]=augmented[pivot],augmented[col]
        divisor=augmented[col][col]
        augmented[col]=[x/divisor for x in augmented[col]]
        for row in range(n):
            if row==col:continue
            factor=augmented[row][col]
            augmented[row]=[x-factor*y for x,y in zip(augmented[row],augmented[col])]
    return [row[-1] for row in augmented]


def project(initial, mass, constraints, target):
    if (len(mass)!=len(initial) or any(x<=0 for x in mass)
            or len(constraints)!=len(target)
            or any(len(row)!=len(initial) for row in constraints)):
        raise ValueError('Invalid constraint dimensions or positive metric')
    schur=[[sum(a*b/m for a,b,m in zip(row,other,mass))
            for other in constraints] for row in constraints]
    residual=[sum(c*u for c,u in zip(row,initial))-value
              for row,value in zip(constraints,target)]
    multiplier=solve(schur,residual)
    velocity=[u-sum(row[i]*value for row,value in zip(constraints,multiplier))/mass[i]
              for i,u in enumerate(initial)]
    return velocity,multiplier


def spread(interpolation, area, traction, volumes):
    """S = Mv^-1 J^T W: force density, with traction in N/m2."""
    return [sum(row[i]*weight*force for row,weight,force in zip(interpolation,area,traction))/volume
            for i,volume in enumerate(volumes)]
