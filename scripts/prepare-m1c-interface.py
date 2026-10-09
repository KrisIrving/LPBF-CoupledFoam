#!/usr/bin/env python3
"""Bounded interface-source matrix on a synthetic pure liquid/vapour pair."""
import importlib.util
import json
import math
import re
from pathlib import Path
import argparse

spec = importlib.util.spec_from_file_location('spatial', Path(__file__).with_name('prepare-m1b-case.py'))
spatial = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spatial)
CASES = ('off', 'prescribed-coarse', 'prescribed-fine', 'prescribed-half-dt', 'kinetic', 'condensation', 'open')


def prepare(case, selection):
    if selection not in CASES:
        raise ValueError('Unknown interface case')
    spatial.prepare(case, 'spatial-evaporation', 'coarse')
    ny = 36 if selection in ('prescribed-fine', 'open') else 18
    dt, end = (5e-10 if selection == 'prescribed-half-dt' else 1e-9), 1e-7
    enabled = selection != 'off'
    model = 'kinetic' if selection in ('kinetic','condensation') else 'prescribed'
    pressure = 120000 if selection == 'condensation' else 80000
    if selection == 'condensation':
        for field in ('p','p_rgh'):
            path=case/'0'/field
            path.write_text(path.read_text().replace('80000','120000'))
    length, floor = 1e-4, 1e-6
    mesh = case/'system/blockMeshDict'
    mesh.write_text(mesh.read_text().replace('(6 6 6)', f'(6 {ny} 6)'))
    # Planar transition; all cell layers are initialized explicitly by setFields.
    # No dependence on blockMesh internal cell ordering.
    layers = []
    for layer in range(ny):
        y = (layer+0.5)*length/ny
        alpha = min(1-floor, max(floor, 0.5*(1-math.tanh((y-length/2)/(length/20)))))
        layers.append(f'''boxToCell {{ box (-1 {layer*length/ny:.17g} -1) (1 {(layer+1)*length/ny:.17g} 1);
fieldValues (volScalarFieldValue alpha.metal1 {alpha:.17g}
volScalarFieldValue alpha.metal1vapour {1-alpha:.17g}); }}''')
    defaults = '\n'.join(f'volScalarFieldValue {p.name} 0' for p in sorted((case/'0').glob('alpha.*')))
    (case/'system/setFieldsDict').write_text(spatial.dictionary('setFieldsDict',
        'defaultFieldValues ('+defaults+');\nregions (\n'+'\n'.join(layers)+'\n);'))
    path = case/'system/fvSolution'
    text = path.read_text().replace('nAlphaSubCycles 2;', 'nAlphaSubCycles 1;')
    text = text.replace('phaseChangeEnabled true;', f'''phaseChangeEnabled {str(enabled).lower()};
interfacePhaseChangeModel {model}; interfaceLiquid metal1; interfaceVapour metal1vapour;
interfaceMassFlux 1; interfaceAccommodation 0.1; interfaceInventoryFraction 0.2;''')
    path.write_text(text)
    spatial.validate_solver_entries(text)
    path = case/'system/controlDict'
    text = path.read_text()
    for name, value in [('endTime',end), ('deltaT',dt), ('maxDeltaT',dt), ('writeInterval',round(end/dt))]:
        text = re.sub(r'\b'+name+r'\s+[^;]+;', f'{name} {value:.17g};', text)
    path.write_text(text+'\ninterfaceDiagnostics true;\n')
    if selection == 'open':
        # y-positive face; reservoir pressure, pressure-driven inlet/outlet U.
        for field in ('p_rgh','U','T')+tuple(p.name for p in (case/'0').glob('alpha.*')):
            path = case/'0'/field
            condition = 'type fixedValue; value uniform 80000;' if field == 'p_rgh' else \
                        'type pressureInletOutletVelocity; value uniform (0 0 0);'
            if field=='T': condition='type inletOutlet; inletValue uniform 4101; value uniform 4101;'
            elif field.startswith('alpha.'):
                inlet=1 if field=='alpha.metal1vapour' else 0
                condition=f'type inletOutlet; inletValue uniform {inlet}; value uniform {inlet};'
            path.write_text(spatial.pair.replace_block(path.read_text(), 'front', condition))
    (case/'m1b-case.json').unlink()
    meta = dict(package='M1C-01', revision=1, selection=selection, model=model, enabled=enabled,
                open_boundary=selection=='open', mesh_cells=[6,ny,6], delta_t_s=dt, end_time_s=end,
                initial_pressure_Pa=pressure, initial_temperature_K=4101,
                expected_planar_area_m2=(1-2*floor)*length**2,
                prescribed_mass_flux_kg_m2_s=1, accommodation=0.1,
                molar_mass_kg_mol=0.05, gas_constant_J_mol_K=8.314,
                reference_pressure_Pa=100000, boiling_temperature_K=4101,
                latent_heat_J_kg=2e6, energy_jump_J_kg=1318067.1651452065,
                comparison_mass_scale_kg=(1-2*floor)*length**2*end,
                thresholds=dict(alpha_bound=1e-6, alpha_sum=1e-6,
                    initial_area_relative=1e-6, max_area_relative=0.05,
                    mass_balance_signal_fraction=0.05, mass_balance_absolute_kg=1e-22,
                    source_pair_relative=1e-12, source_energy_relative=1e-12,
                    limiter_fraction=1e-4, required_pressure_rise_Pa=1,
                    minimum_speed_m_s=1e-5),
                scope='Interface source code verification, pure vapour and equal Cv synthetic fixture. '
                      'No Knudsen exit momentum, shielding gas, resolved particles, general energy closure or material validation.')
    (case/'m1c-case.json').write_text(json.dumps(meta, indent=2)+'\n')
    return meta


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',type=Path)
    parser.add_argument('selection',choices=CASES)
    args=parser.parse_args()
    prepare(args.case,args.selection)
