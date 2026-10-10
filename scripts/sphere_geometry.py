"""Independent Cartesian sphere geometry reference; not a CFD discretization.

Circle/rectangle intersections are analytic. Sphere/box volumes integrate those
areas with adaptive Simpson quadrature, split at geometric transition points.
All integration is normalized by the real radius; no radius/volume fitting.
Only Python's standard library is required. Units follow the input coordinates.
"""
import math


def disk_rectangle(radius, lo, hi):
    """Area of a centred disk intersecting a 2D rectangle."""
    if radius <= 0:
        return 0.0

    def primitive(x):
        x = min(radius, max(0.0, x))
        return .5*(x*math.sqrt(max(0.0, radius*radius-x*x))
                   + radius*radius*math.asin(x/radius))

    def quadrant(x, y):
        a, b = min(abs(x), radius), min(abs(y), radius)
        split = min(a, math.sqrt(max(0.0, radius*radius-b*b)))
        value = b*split + primitive(a)-primitive(split)
        return math.copysign(1.0, x)*math.copysign(1.0, y)*value

    area = (quadrant(hi[0], hi[1])-quadrant(lo[0], hi[1])
            -quadrant(hi[0], lo[1])+quadrant(lo[0], lo[1]))
    return min((hi[0]-lo[0])*(hi[1]-lo[1]), max(0.0, area))


def _integrate(function, a, b, tolerance, depth=22):
    def simpson(a, b, fa, fm, fb):
        return (b-a)*(fa+4*fm+fb)/6

    def refine(a, b, fa, fm, fb, whole, eps, remaining):
        mid = (a+b)/2
        fl, fr = function((a+mid)/2), function((mid+b)/2)
        left = simpson(a, mid, fa, fl, fm)
        right = simpson(mid, b, fm, fr, fb)
        correction = (left+right-whole)/15
        if abs(correction) <= eps:
            return left+right+correction, abs(correction)
        if remaining == 0:
            raise ArithmeticError('Sphere-volume quadrature exhausted its depth')
        lv, le = refine(a, mid, fa, fl, fm, left, eps/2, remaining-1)
        rv, re = refine(mid, b, fm, fr, fb, right, eps/2, remaining-1)
        return lv+rv, le+re

    fa, fm, fb = function(a), function((a+b)/2), function(b)
    return refine(a, b, fa, fm, fb, simpson(a,b,fa,fm,fb), tolerance, depth)


def sphere_box(centre, radius, lo, hi, fraction_tolerance=1e-8):
    """Solid volume, six solid face areas, integrated outward surface normal.

    Face order: x-,x+,y-,y+,z-,z+. Fluid aperture is 1-solid_face_fraction.
    normal_area is integral(n dA) on the actual curved surface inside the cell,
    obtained by closure with planar solid faces. It is NOT scalar surface area.
    Quadrature error is an estimate, not a rigorous bound or CFD error estimate.
    """
    if (radius <= 0 or len(centre) != 3 or len(lo) != 3 or len(hi) != 3
            or not 0 < fraction_tolerance <= 1e-3
            or not all(math.isfinite(x) for x in (*centre, radius, *lo, *hi))
            or any(b <= a for a,b in zip(lo,hi))):
        raise ValueError('Finite positive sphere and nondegenerate box required')
    lower = [(x-c)/radius for x,c in zip(lo,centre)]
    upper = [(x-c)/radius for x,c in zip(hi,centre)]
    box_volume = math.prod(b-a for a,b in zip(lo,hi))
    minimum2 = sum((a if a > 0 else b if b < 0 else 0)**2
                   for a,b in zip(lower,upper))
    maximum2 = sum(max(a*a,b*b) for a,b in zip(lower,upper))
    face_areas = []
    for axis in range(3):
        other = [j for j in range(3) if j != axis]
        for position in (lower[axis], upper[axis]):
            face_areas.append(radius**2*disk_rectangle(
                math.sqrt(max(0.0,1-position*position)),
                [lower[j] for j in other], [upper[j] for j in other]))
    error = 0.0
    if minimum2 >= 1:
        volume, classification = 0.0, 'outside'
    elif maximum2 <= 1:
        volume, classification = box_volume, 'inside'
    else:
        classification = 'cut'
        a, b = max(-1.0,lower[0]), min(1.0,upper[0])
        transitions = {a,b}
        # Circle reaches each edge/corner of the transverse rectangle.
        for y in (0.0,lower[1],upper[1]):
            for z in (0.0,lower[2],upper[2]):
                distance2 = y*y+z*z
                if distance2 < 1:
                    point = math.sqrt(1-distance2)
                    transitions.update(x for x in (-point,point) if a < x < b)
        knots = sorted(transitions)
        eps = fraction_tolerance*box_volume/radius**3/(len(knots)-1)
        def section(x):
            return disk_rectangle(math.sqrt(max(0.0,1-x*x)),lower[1:],upper[1:])
        pieces = [_integrate(section,x,y,eps) for x,y in zip(knots,knots[1:])]
        volume = radius**3*sum(v for v,e in pieces)
        error = radius**3*sum(e for v,e in pieces)
        volume = min(box_volume,max(0.0,volume))
    fractions = []
    for axis in range(3):
        face_size = math.prod(hi[j]-lo[j] for j in range(3) if j != axis)
        fractions.extend(face_areas[2*axis+k]/face_size for k in (0,1))
    return {'classification': classification, 'volume': volume,
            'solid_fraction': volume/box_volume, 'solid_face_areas': face_areas,
            'fluid_face_apertures': [1-x for x in fractions],
            'normal_area': [face_areas[2*j]-face_areas[2*j+1] for j in range(3)],
            'estimated_volume_error': error}


def segment_sphere(start, end, centre, radius):
    """Exact straight segment intersections and outward unit normals.

    Tangency is returned once; endpoints included. No crossing implies [] even
    for a segment fully inside. This routine does not classify cell interiors.
    """
    if (radius <= 0 or any(len(x) != 3 for x in (start,end,centre))
            or not all(math.isfinite(x) for x in (*start,*end,*centre,radius))):
        raise ValueError('Finite positive sphere required')
    d = [b-a for a,b in zip(start,end)]
    r = [a-c for a,c in zip(start,centre)]
    a = sum(x*x for x in d)
    if a == 0:
        raise ValueError('Nonzero segment required')
    b, c = 2*sum(x*y for x,y in zip(d,r)), sum(x*x for x in r)-radius*radius
    discriminant = b*b-4*a*c
    if discriminant < 0:
        return []
    # Stable roots avoid cancellation for a crossing close to an endpoint.
    q = -.5*(b+math.copysign(math.sqrt(discriminant),b))
    roots = {-b/(2*a)} if discriminant == 0 else {q/a,c/q}
    result = []
    for t in sorted(roots):
        if 0 <= t <= 1:
            point = [x+t*y for x,y in zip(start,d)]
            result.append({'t':t, 'point':point,
                           'normal':[(x-y)/radius for x,y in zip(point,centre)]})
    return result
