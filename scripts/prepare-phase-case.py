#!/usr/bin/env python3
"""Prepare a homogeneous closed-box synthetic pair from copied Test1."""
import argparse
import json
import re
from pathlib import Path

PATCHES = ('top', 'bottom', 'right', 'left', 'front', 'back')


def field(path, value, condition):
    text = path.read_text()
    text, count = re.subn(r'internalField\s+uniform\s+[^;]+;',
                         f'internalField uniform {value};', text)
    if count != 1:
        raise ValueError(f'Expected uniform initial field: {path}')
    text = text[:text.index('boundaryField')] + 'boundaryField\n{\n'
    for patch in PATCHES:
        text += f'    {patch}\n    {{\n        {condition}\n    }}\n'
    path.write_text(text+'}\n')


def prepare(case, direction):
    thermo = (case/'constant/thermophysicalProperties').read_text()
    tables = {}
    for name in ('boils', 'LatentHeatGas'):
        table = re.search(r'\b'+name+r'\s*\((.*?)\);', thermo, re.S).group(1)
        tables[name] = float(re.search(r'\(metal1vapour\s+metal1\)\s+([\deE.+-]+)', table).group(1))
    p0 = float(re.search(r'\bP0\s+([\deE.+-]+)\s*;', thermo).group(1))
    if (tables['boils'], tables['LatentHeatGas'], p0) != (4101, 10, 100000):
        raise ValueError('Synthetic test requires the documented upstream pair parameters')
    # At T=4101 K, this upstream pair has Psat=P0=100000 Pa.
    pressure = {'evaporation': 80000, 'condensation': 120000, 'equilibrium': 100000}[direction]
    for path in (case/'0').glob('alpha.*'):
        value = '0.5' if path.name in ('alpha.metal1', 'alpha.metal1vapour') else '0'
        field(path, value, 'type zeroGradient;')
    field(case/'0/T', '4101', 'type zeroGradient;')
    field(case/'0/U', '(0 0 0)', 'type noSlip;')
    field(case/'0/p', str(pressure), f'type calculated; value uniform {pressure};')
    field(case/'0/p_rgh', str(pressure), f'type fixedFluxPressure; value uniform {pressure};')
    p = case/'constant/g'
    text, count = re.subn(r'value\s+\([^)]*\)\s*;', 'value (0 0 0);', p.read_text())
    if count != 1:
        raise ValueError('Expected gravity vector')
    p.write_text(text)
    (case/'constant/timeVsLaserPower').write_text('(\n (0 0)\n (1 0)\n)\n')
    (case/'phase-case.json').write_text(json.dumps(dict(
        direction=direction, temperature_K=4101, initial_pressure_Pa=pressure,
        initial_saturation_pressure_Pa=100000, latent_heat_J_per_kg=10,
        liquid_alpha=0.5, vapour_alpha=0.5,
        interpretation='Synthetic homogeneous closure test; upstream EOS and latent parameters retained. '
                       'Remaining phases are inactive at initialization, not removed. '
                       'Initial driving direction may change as pressure and temperature evolve.'
    ), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    parser.add_argument('direction', choices=('evaporation', 'condensation', 'equilibrium'))
    args = parser.parse_args()
    prepare(args.case, args.direction)
