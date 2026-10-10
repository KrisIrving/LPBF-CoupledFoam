"""Generate the bounded user-side M2A-01 mechanical verification matrix."""
import argparse
import json
from pathlib import Path

CASES = {'fixed-coarse': ('fixed', 36), 'fixed-fine': ('fixed', 48),
         'rotate-coarse': ('rotate', 36), 'rotate-fine': ('rotate', 48),
         'translate': ('translate', 36), 'free': ('free', 36),
         'free-restart': ('free', 36), 'free-mpi2': ('free', 36)}


def foam_header(kind, name):
    return f'FoamFile {{ version 2.0; format ascii; class {kind}; object {name}; }}\n'


def write_dem(root, mode, restart=False):
    lines = ['units si', 'atom_style sphere', 'atom_modify map array', 'newton off',
             'boundary f f f', 'communicate single vel yes', 'processors * 1 1',
             'read_restart 0.002/dem.restart' if restart else 'read_data particle.data',
             'fix Y all property/global youngsModulus peratomtype 1e7',
             'fix P all property/global poissonsRatio peratomtype 0.3',
             'fix R all property/global coefficientRestitution peratomtypepair 1 0.5',
             'fix F all property/global coefficientFriction peratomtypepair 1 0',
             'pair_style gran model hertz tangential history', 'pair_coeff * *',
             'neighbor 0.002 bin', 'neigh_modify delay 0 every 1 check yes',
             'fix integrate all nve/sphere', 'timestep 0.00001', 'thermo 100',
             'thermo_style custom step atoms time ke', 'thermo_modify lost error',
             'variable demoTime equal step*dt', 'run 0']
    (root / 'input.dem').write_text('\n'.join(lines)+'\n')


def control(root, start=0., end=0.004):
    (root / 'system/controlDict').write_text(foam_header('dictionary', 'controlDict') + f'''
application resolvedParticleFoam;
startFrom startTime;
startTime {start:.12g};
stopAt endTime;
endTime {end:.12g};
deltaT 0.0001;
writeControl timeStep;
writeInterval 20;
writeFormat ascii;
writePrecision 17;
timeFormat general;
timePrecision 12;
runTimeModifiable false;
''')


def prepare(root, selection):
    mode, n = CASES[selection]
    for d in ('system', 'constant', '0'):
        (root/d).mkdir(parents=True, exist_ok=True)
    control(root)
    (root / 'system/blockMeshDict').write_text(foam_header('dictionary', 'blockMeshDict')+f'''
scale 1;
vertices ((-0.06 -0.06 -0.06) (0.06 -0.06 -0.06) (0.06 0.06 -0.06) (-0.06 0.06 -0.06)
          (-0.06 -0.06 0.06) (0.06 -0.06 0.06) (0.06 0.06 0.06) (-0.06 0.06 0.06));
blocks (hex (0 1 2 3 4 5 6 7) ({n} {n} {n}) simpleGrading (1 1 1));
edges ();
boundary (outer {{ type patch; faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }});
mergePatchPairs ();
''')
    (root / 'system/fvSchemes').write_text(foam_header('dictionary', 'fvSchemes')+'''
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(advecting,U) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; }
''')
    (root / 'system/fvSolution').write_text(foam_header('dictionary', 'fvSolution')+'''
solvers
{
    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0; }
    p { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }
}
''')
    (root / 'system/decomposeParDict').write_text(foam_header('dictionary', 'decomposeParDict')+'''
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
''')
    (root / 'constant/mechanicalProperties').write_text(foam_header('dictionary', 'mechanicalProperties')+f'''
mode {mode};
rho 1;
nu 0.1;
demDeltaT 0.00001;
penalty 1000;
quadrature 4;
correctors 8;
maxCorrectors 96;
momentumImpulseTolerance 1e-11;
''')
    (root / '0/U').write_text(foam_header('volVectorField', 'U')+'''
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField { outer { type fixedValue; value uniform (0 0 0); } }
''')
    (root / '0/p').write_text(foam_header('volScalarField', 'p')+'''
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField { outer { type fixedFluxPressure; value uniform 0; } }
''')
    speed = 0.3 if mode == 'translate' else 0.01 if mode == 'free' else 0.
    centre = -0.0005 if mode in {'translate', 'free'} else 0.
    omega = 0.1 if mode == 'rotate' else 0.
    (root / 'particle.data').write_text(f'''M2A noncontact sphere; intentionally sparse global ID

1 atoms
1 atom types

-0.06 0.06 xlo xhi
-0.06 0.06 ylo yhi
-0.06 0.06 zlo zhi

Atoms

7 1 0.02 1000 {centre} 0 0

Velocities

7 {speed} 0 0 0 0 {omega}
''')
    write_dem(root, mode)
    metadata = {'package': 'M2A-01', 'selection': selection, 'mode': mode, 'mesh_n': n,
                'rho': 1., 'nu': 0.1, 'particle_density': 1000., 'radius': 0.01, 'id': 7,
                'initial_centre': [centre, 0., 0.], 'initial_velocity': [speed, 0., 0.],
                'initial_omega': [0., 0., omega], 'delta_t': 1e-4, 'dem_delta_t': 1e-5,
                'end_time': 0.004, 'expected_rows': 40,
                'reference_flow_speed': 0.01, 'reference_rotation_speed': 0.1,
                'solver_controls': {'min_correctors': 8, 'max_correctors': 96,
                                    'momentum_impulse_L1_kg_m_s': 1e-11},
                'thresholds': {'volume_relative': 0.05, 'slip_relative': 0.05,
                               'divergence_per_s': 1e-6, 'momentum_residual_kg_m_s': 1e-10,
                               'angular_impulse_residual_kg_m2_s': 1e-12,
                               'clock_s': 1e-12, 'stokes_force_relative': 0.35,
                               'stokes_torque_relative': 0.35,
                               'force_decomposition_N': 1e-12, 'torque_decomposition_N_m': 1e-14,
                               'refinement_error_increase': 0.02,
                               'comparison_position_m': 1e-8, 'comparison_velocity_m_s': 1e-8,
                               'comparison_rotation_per_s': 1e-8, 'comparison_force_N': 1e-8,
                               'comparison_torque_N_m': 1e-10, 'comparison_fluid_momentum_kg_m_s': 1e-10},
                'scope': 'Experimental single noncontact sphere, cubic mesh, implicit volume '
                         'constraint, frozen geometry per window and virtual-fluid inertia audit. '
                         'No laser, evaporation, heat, collisions or production scaling.'}
    (root / 'M2A_META.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', type=Path)
    parser.add_argument('--selection', choices=CASES)
    parser.add_argument('--restart', action='store_true')
    parser.add_argument('--half', action='store_true')
    args = parser.parse_args()
    if args.restart:
        meta = json.loads((args.case/'M2A_META.json').read_text())
        control(args.case, start=0.002)
        write_dem(args.case, meta['mode'], restart=True)
    elif args.half:
        control(args.case, end=0.002)
    elif args.selection:
        prepare(args.case, args.selection)
    else:
        parser.error('Provide selection, half or restart')
