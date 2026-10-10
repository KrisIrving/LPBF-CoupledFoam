"""Six-fixture geometry reference audit; no OpenFOAM/LIGGGHTS required.

This does not alter the mechanical solver or replace its failed force gate.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import time

from sphere_geometry import sphere_box

THRESHOLDS = {'volume_relative_error':1e-6, 'normal_closure_over_area':1e-10,
              'shared_face_relative_difference':1e-12}


def fixture(n, phase):
    start = time.perf_counter()
    h, radius, origin = .12/n,.01,-.06
    centre = [phase*h]*3
    exact = 4*math.pi*radius**3/3
    ranges = [range(math.floor((c-radius-origin)/h),
                    math.ceil((c+radius-origin)/h)) for c in centre]
    cells, volume, estimate, sample_volume = {},0.,0.,0.
    support_reference = support_sample = 0
    normal = [0.,0.,0.]
    for index in itertools.product(*ranges):
        lo = [origin+i*h for i in index]
        hi = [origin+(i+1)*h for i in index]
        geometry = sphere_box(centre,radius,lo,hi)
        cells[index] = geometry
        volume += geometry['volume']
        estimate += geometry['estimated_volume_error']
        if geometry['volume'] > h**3*1e-12:
            support_reference += 1
        normal = [a+b for a,b in zip(normal,geometry['normal_area'])]
        hits = 0
        for sample in itertools.product(range(4),repeat=3):
            point = [x+(s+.5)*h/4 for x,s in zip(lo,sample)]
            hits += sum((x-c)**2 for x,c in zip(point,centre)) <= radius**2
        sample_volume += h**3*hits/64
        support_sample += hits > 0
    shared = 0.
    for index, cell in cells.items():
        for axis in range(3):
            neighbour = list(index); neighbour[axis] += 1
            if tuple(neighbour) in cells:
                other = cells[tuple(neighbour)]
                shared = max(shared,abs(cell['solid_face_areas'][2*axis+1]
                                       -other['solid_face_areas'][2*axis])/h**2)
    metrics = {'volume_relative_error':abs(volume/exact-1),
               'normal_closure_over_area':math.sqrt(sum(x*x for x in normal))/(4*math.pi*radius**2),
               'shared_face_relative_difference':shared,
               'estimated_volume_error_relative':estimate/exact,
               'old_q4_volume_relative_error':abs(sample_volume/exact-1),
               'old_q4_support_volume_ratio':support_sample*h**3/exact,
               'intersection_cell_support_volume_ratio':support_reference*h**3/exact,
               'elapsed_seconds':time.perf_counter()-start}
    checks = {key:metrics[key] <= limit for key,limit in THRESHOLDS.items()}
    return {'n':n,'diameter_over_h':2*radius/h,'centre_phase_h':phase,
            'metrics':metrics,'checks':checks,'passed':all(checks.values())}


def run():
    cases = {f'n{n}-phase{phase}':fixture(n,phase)
             for n in (36,48,64) for phase in (0.,.25)}
    return {'package':'M2A-02-geometry-reference','thresholds':THRESHOLDS,
            'cases':cases,'passed':all(x['passed'] for x in cases.values()),
            'scope':'Offline geometry reference only; no new CFD force, velocity '
                    'constraint, moving GCL, MPI mapping or LPBF validation.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for name,case in result['cases'].items():
        print(name,'PASS' if case['passed'] else 'FAIL',case['metrics'])
    raise SystemExit(0 if result['passed'] else 1)
