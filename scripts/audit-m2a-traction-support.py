"""Geometry-only audit of complete-Q source support and traction probe stencils.

An overlap is potential source support, not evidence of actual nonzero forcing.
No flow, pressure projection or physical acceptance is computed.
"""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import surface_reference as ref

spec = importlib.util.spec_from_file_location('sampling', Path(__file__).with_name('audit-m2a-complete-wall.py'))
sampling = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sampling)


def overlap(row, support):
    return sum(weight for cell, weight in row.items() if cell in support)


def fixture(n, phase):
    h = .12/n
    radius = .01
    origin = [-.06+h/2]*3
    centre = [phase*h]*3
    quadrature = list(ref.sphere_quadrature(radius))
    wall_rows = [sampling.sampling_row([c+radius*x for c, x in zip(centre, normal)], origin, h)
                 for normal, _ in quadrature]
    support = {cell for row in wall_rows for cell, weight in row.items() if weight != 0}
    radial = [math.sqrt(sum((o+i*h-c)**2 for o, i, c in zip(origin, cell, centre))) for cell in support]
    probes = {}
    for distance in (2, 3, 4):
        fractions = []
        for normal, area in quadrature:
            point = [c+(radius+distance*h)*x for c, x in zip(centre, normal)]
            fractions.append((overlap(sampling.sampling_row(point, origin, h), support), area))
        probes[str(distance)] = {
            'nodes_with_support_overlap': sum(value > 0 for value, _ in fractions),
            'max_interpolation_weight_in_support': max(value for value, _ in fractions),
            'area_mean_interpolation_weight_in_support': sum(value*area for value, area in fractions)/sum(area for _, area in fractions)}
    return {'mesh_n': n, 'phase_h': phase, 'support_cells': len(support),
            'exterior_support_cells': sum(r > radius for r in radial),
            'maximum_support_distance_outside_wall_h': (max(radial)-radius)/h,
            'traction_probes_distance_h': probes}


def run():
    return {'scope': __doc__, 'cases': [fixture(n, phase) for n in (36, 48, 64) for phase in (0., .25)],
            'physical_pass': 'not_assessed',
            'limitation': 'Stencil union bounds possible forcing only; extrapolation consistency on the numerical field remains unverified.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False)+'\n')
