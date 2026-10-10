"""Offline formula oracle, not the C++ runtime or a CFD simulation."""
import itertools
import math


def trilinear(point,origin,h,values):
    coordinate=[(x-o)/h for x,o in zip(point,origin)]
    lower=[math.floor(x) for x in coordinate]
    fraction=[x-i for x,i in zip(coordinate,lower)]
    result=[0.]*len(next(iter(values.values())))
    for offsets in itertools.product((0,1),repeat=3):
        weight=math.prod(f if i else 1-f for f,i in zip(fraction,offsets))
        value=values[tuple(i+j for i,j in zip(lower,offsets))]
        result=[a+weight*b for a,b in zip(result,value)]
    return result


def ghost_target(wall,image,distance,h):
    return [a+distance/(2*h)*(b-a) for a,b in zip(wall,image)]


def radial_derivative(wall,a,b,c,h):
    return [(-13*x/12+3*y-8*z/3+3*t/4)/h for x,y,z,t in zip(wall,a,b,c)]


def quadratic_ghost_target(wall,a,b,distance,h):
    t=distance/h
    return [(t-2)*(t-3)/6*x-t*(t-3)/2*y+t*(t-2)/3*z for x,y,z in zip(wall,a,b)]


def pressure_wall(a,b,c):
    return 6*a-8*b+3*c


def sphere_quadrature(radius,nz=12,nphi=24):
    for iz in range(nz):
        z=math.cos(math.pi*(iz+.75)/(nz+.5))
        for iteration in range(32):
            a,b=1.,0.
            for degree in range(1,nz+1):
                a,b=((2*degree-1)*z*a-(degree-1)*b)/degree,a
            derivative=nz*(z*a-b)/(z*z-1)
            change=a/derivative;z-=change
            if abs(change)<1e-15:
                break
        else:
            raise ArithmeticError('Quadrature root budget exhausted')
        weight=2/((1-z*z)*derivative**2)*2*math.pi/nphi*radius**2
        for ip in range(nphi):
            phi=2*math.pi*(ip+.5)/nphi
            yield [math.sqrt(1-z*z)*math.cos(phi),math.sqrt(1-z*z)*math.sin(phi),z],weight


def cross(a,b):
    return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]


def traction(normal,radial,omega,pressure,rho,nu):
    spin=cross(omega,normal)
    projection=sum(a*b for a,b in zip(normal,radial))
    return ([-rho*pressure*x for x in normal],
            [rho*nu*(a+n*projection-s) for a,n,s in zip(radial,normal,spin)])
