"""Offline audit of analytic Dirichlet face-centre flux; no CFD execution.

A closed box with prescribed velocity needs sum(Uf.Sf)=0, independently of
pressure tolerance. Compare face-centre integration with 2x2 Gauss integration.
The reference-cell estimate is a diagnostic, not a pressure-solver simulation.
"""
import argparse
import itertools
import json
import math
from pathlib import Path


def velocity(point, centre, mode):
    r = [x-c for x, c in zip(point, centre)]
    distance = math.sqrt(sum(x*x for x in r))
    a = .01/distance
    if mode == 'rotate':
        return [-.1*r[1]*a**3, .1*r[0]*a**3, 0.]
    normal = [x/distance for x in r]
    return [(1-.75*a-.25*a**3)*u + (-.75*a+.75*a**3)*.01*normal[0]*v
            for u, v in zip((.01, 0., 0.), normal)]


def fixture(n, phase, mode, gauss=False):
    h = .12/n
    centre = [phase*h]*3
    offsets = (-1/math.sqrt(3), 1/math.sqrt(3)) if gauss else (0.,)
    contributions = []
    for axis, side in itertools.product(range(3), (-1, 1)):
        tangent = [i for i in range(3) if i != axis]
        for i, j in itertools.product(range(n), repeat=2):
            values = []
            for s, t in itertools.product(offsets, repeat=2):
                point = [0.]*3
                point[axis] = side*.06
                point[tangent[0]] = -.06+(i+.5)*h+s*h/2
                point[tangent[1]] = -.06+(j+.5)*h+t*h/2
                values.append(side*velocity(point, centre, mode)[axis])
            contributions.append(math.fsum(values)*h*h/len(values))
    flux = math.fsum(contributions)
    return {'mesh_n': n, 'phase_h': phase, 'mode': mode,
            'quadrature': 'gauss2x2' if gauss else 'faceCentre',
            'net_boundary_flux_m3_s': flux,
            'reference_cell_divergence_estimate_per_s': abs(flux)/h**3,
            'box_mean_divergence_lower_bound_per_s': abs(flux)/.12**3}


def run():
    return {'kind': 'analytic_boundary_quadrature_only',
            'limitations': ['No pressure matrix or CFD execution',
                            'Gauss integration is a diagnostic, not an applied boundary fix'],
            'cases': [fixture(n, phase, mode, gauss)
                      for mode in ('fixed', 'rotate') for n in (36, 48, 64)
                      for phase in (0., .25) for gauss in (False, True)]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False)+'\n')
