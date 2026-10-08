#!/usr/bin/env python3
"""Generate spatial synthetic integration cases; no material reproduction claim."""
import argparse
import importlib.util
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('pair', ROOT/'scripts/prepare-phase-case.py')
pair = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pair)
FAMILIES = ('advection', 'spatial-evaporation', 'laser-heat', 'laser-phase')


def dictionary(name, body):
    return f'FoamFile\n{{\n version 2.0; format ascii; class dictionary; object {name};\n}}\n'+body+'\n'


def prepare(case, family, level, profile='short'):
    if family not in FAMILIES or level not in ('coarse', 'fine') or profile not in ('short', 'verification'):
        raise ValueError('Unsupported integrated test selection')
    case.mkdir(parents=True, exist_ok=False)
    for part in ('constant', 'system'):
        shutil.copytree(ROOT/'tutorials/compressiblelaserbeamFoam/Test1'/part, case/part)
    shutil.copytree(ROOT/'tutorials/compressiblelaserbeamFoam/Test1/initial', case/'0')
    pair.prepare(case, 'evaporation' if family == 'spatial-evaporation' else 'equilibrium', 'consistent')
    (case/'phase-case.json').unlink()  # Uniform-case audit must not be applied here.
    advection = family == 'advection'
    laser_on = family.startswith('laser-')
    active = family in ('spatial-evaporation', 'laser-phase')
    nx, ny, nz = ((10 if level == 'coarse' else 20), 2, 2) if advection else ((6,)*3 if level == 'coarse' else (12,)*3)
    length = 1e-4
    pressure = 80000 if family == 'spatial-evaporation' else 100000
    dt = 1e-6 if advection else 5e-10
    end = (2e-5 if advection else 1e-8)*(2 if profile == 'verification' else 1)
    velocity = 1 if advection else 0
    power = 1 if laser_on else 0
    offset = 0 if advection else 1318067.1651452065
    if advection:
        # Matched rho/Cv removes density-contrast physics from the transport reference.
        path = case/'constant/thermophysicalProperties.metal1vapour'
        text = re.sub(r'\bequationOfState\s+perfectGas;', 'equationOfState rhoConst;', path.read_text())
        path.write_text(pair.replace_block(text, 'equationOfState', 'rho 5000;'))
    path = case/'constant/thermophysicalProperties'
    text = path.read_text()
    for table in ('sigmas', 'dsigmadT', 'interfaceDiffusion')+ (('interfaceCompression',) if advection else ()):
        match = re.search(r'\b'+table+r'\s*\((.*?)\);', text, re.S)
        values = re.sub(r'(\([^()]+\)\s+)[\deE.+-]+', r'\g<1>0', match.group(1))
        text = text[:match.start(1)]+values+text[match.end(1):]
    text = re.sub(r'(\bmetal1vapour\s+)1318067[^;]*;', r'\g<1>'+str(offset)+';', text)
    path.write_text(text)
    mesh = (case/'system/blockMeshDict').read_text().replace('scale 0.001;', 'scale 0.00001;')
    mesh = re.sub(r'\(5 5 5\)', f'({nx} {ny} {nz})', mesh)
    if advection:
        for patch, neighbour in (('left', 'right'), ('right', 'left')):
            mesh = re.sub(r'(\b'+patch+r'\s*\{\s*)type patch;',
                          r'\g<1>type cyclic; neighbourPatch '+neighbour+';', mesh)
    (case/'system/blockMeshDict').write_text(mesh)
    for path in (case/'0').iterdir():
        if path.name not in ('T', 'U', 'p', 'p_rgh', 'Laser_boundary') and not path.name.startswith('alpha.'):
            continue
        name = path.name
        value = '(1 0 0)' if name == 'U' and advection else '(0 0 0)' if name == 'U' else '4101' if name == 'T' else str(pressure) if name in ('p', 'p_rgh') else '0'
        condition = 'type slip;' if name == 'U' and advection else 'type noSlip;' if name == 'U' else 'type calculated; value uniform '+str(pressure)+';' if name == 'p' else 'type fixedFluxPressure; value uniform '+str(pressure)+';' if name == 'p_rgh' else 'type zeroGradient;'
        pair.field(path, value, condition)
        content = path.read_text()
        if advection:
            for patch in ('left', 'right'):
                content = pair.replace_block(content, patch, 'type cyclic;')
            if name == 'p_rgh':
                # This upstream pressure branch has no all-incompressible reference
                # handling. Fix one reference plane; uniform translation has zero flux there.
                content = pair.replace_block(content, 'top', 'type fixedValue; value uniform 100000;')
        if name == 'Laser_boundary':
            content = re.sub(r'object\s+T;', 'object Laser_boundary;', content)
            for patch in pair.PATCHES:
                if advection and patch in ('left', 'right'):
                    continue
                content = pair.replace_block(content, patch, 'type fixedValue; value uniform '+('1' if patch == 'back' else '-1')+';')
        path.write_text(content)
    slab = '(2e-5 -1 -1) (5e-5 1 1)' if advection else '(-1 5e-5 -1) (1 1 1)'
    defaults = '\n'.join('volScalarFieldValue '+p.name+' '+('1' if p.name == 'alpha.metal1vapour' else '0') for p in sorted((case/'0').glob('alpha.*')))
    (case/'system/setFieldsDict').write_text(dictionary('setFieldsDict', f'''
defaultFieldValues ({defaults});
regions (boxToCell {{ box {slab}; fieldValues (
volScalarFieldValue alpha.metal1 1
volScalarFieldValue alpha.metal1vapour 0
); }});
'''))
    reference = f'''spatialReference {{ cellWidthX {length/nx:.17g}; velocityX 1;
slabMinX 2e-5; slabMaxX 5e-5; }}''' if advection else ''
    (case/'system/controlDict').write_text(dictionary('controlDict', f'''
application compressibleLaserbeamFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {end:.17g};
deltaT {dt:.17g}; adjustTimeStep false; maxCo 0.5; maxAlphaCo 0.5; maxDeltaT {dt:.17g};
writeControl timeStep; writeInterval {round(end/dt)}; purgeWrite 0;
writeFormat ascii; writePrecision 15; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false;
continuumDiagnostics true; couplingDiagnostics true;
temperatureBudgetDiagnostics true; spatialDiagnostics true;
{reference}
'''))
    (case/'system/fvSolution').write_text(dictionary('fvSolution', f'''
solvers
{{
 "alpha.*" {{ nAlphaSubCycles {1 if advection else 2}; cAlpha {0 if advection else 1};
 averagePhaseChangeSources true; protectAllPhaseOldTimes true;
 commonLatentHeatSource true; referenceInternalEnergySource true;
 phaseChangeEnabled {str(active).lower()}; }}
 "rho.*" {{ solver diagonal; }}
 p_rgh {{ solver GAMG; tolerance 1e-8; relTol 0; smoother GaussSeidel; maxIter 1000; }}
 p_rghFinal {{ $p_rgh; }}
 T {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; maxIter 1000; }}
 TFinal {{ $T; }}
 U {{ solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0; maxIter 1000; }}
 UFinal {{ $U; }}
}}
MELTING {{ minTempCorrector 1; maxTempCorrector 20; epsilonTolerance 1e-6;
 epsilonRelaxation 0.9; damperSwitch false; }}
PIMPLE {{ nOuterCorrectors 5; nCorrectors 3; nNonOrthogonalCorrectors 0; }}
relaxationFactors {{ equations {{ "U.*" 1; }} }}
'''))
    path = case/'constant/dynamicMeshDict'
    path.write_text(re.sub(r'dynamicFvMesh\s+\w+;', 'dynamicFvMesh staticFvMesh;', path.read_text()))
    path = case/'constant/LaserProperties'
    path.write_text(re.sub(r'laserRadius\s+[^;]+;', 'laserRadius 2.5e-5;', path.read_text()))
    (case/'constant/timeVsLaserPosition').write_text('(\n(0 (5e-5 0 5e-5))\n(1 (5e-5 0 5e-5))\n)\n')
    (case/'constant/timeVsLaserPower').write_text(f'(\n(0 {power})\n(1 {power})\n)\n')
    meta = dict(family=family, level=level, profile=profile, cells=nx*ny*nz,
                mesh_cells=[nx,ny,nz], length_m=length, initial_temperature_K=4101,
                initial_pressure_Pa=pressure, delta_t_s=dt, end_time_s=end,
                phase_change_enabled=active, laser_power_w=power, velocity_x_m_s=velocity,
                vapour_reference_energy_offset_J_per_kg=offset,
                thresholds=dict(alpha_bound_tolerance=1e-6, alpha_sum_tolerance=1e-6,
                    closed_mass_relative_tolerance=1e-5, advection_mean_abs_error=0.12,
                    advection_phase_volume_relative_tolerance=1e-5,
                    minimum_active_temperature_source_kg_k_per_s=1e-16),
                scope='Spatial synthetic integration tests only. Equal Cv; matched density for '
                      'advection, rhoConst/perfectGas otherwise. No real material, melting or '
                      'literature validation. Native-energy change and absorbed laser energy are '
                      'observations, not a universal energy conservation gate.')
    (case/'m1b-case.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',type=Path)
    parser.add_argument('family',choices=FAMILIES)
    parser.add_argument('level',choices=('coarse','fine'))
    parser.add_argument('--profile',choices=('short','verification'),default='short')
    args=parser.parse_args()
    prepare(args.case,args.family,args.level,args.profile)
