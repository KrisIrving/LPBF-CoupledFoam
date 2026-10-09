"""Offline feedback-contract tests; no OpenFOAM execution."""
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


def module(name,file):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/file)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


prepare=module('prepare','prepare-m1c-interface.py')
summary=module('summary','summarize-m1c-interface.py')


def record(prefix,**values):
    return prefix+' '+' '.join(f'{k}={v}' for k,v in values.items())+'\n'


def feedback(case,selection='prescribed-coarse'):
    meta=prepare.prepare(case,selection)
    meta['end_time_s']=2*meta['delta_t_s']
    meta['comparison_mass_scale_kg']=meta['expected_planar_area_m2']*meta['end_time_s']
    (case/'m1c-case.json').write_text(json.dumps(meta))
    names=[p.name.removeprefix('alpha.') for p in sorted((case/'0').glob('alpha.*'))]
    rate=meta['expected_planar_area_m2'] if meta['enabled'] else 0
    if meta['model']=='kinetic':
        rate*=meta['accommodation']*(meta['reference_pressure_Pa']-meta['initial_pressure_Pa'])*math.sqrt(
            meta['molar_mass_kg_mol']/(2*math.pi*meta['gas_constant_J_mol_K']*meta['initial_temperature_K']))
    log=''
    for index in range(3):
        time=f'{index*meta["delta_t_s"]:.15g}'
        transfer=rate*float(time)
        # Two trial sources per step; only the final trial is integrated.
        if index:
            for trial in (rate if index==1 else 2*rate,rate):
                log+=record('M1C_SOURCE',time=time,dt=meta['delta_t_s'],model=meta['model'],enabled=str(meta['enabled']).lower(),
                            areaM2=meta['expected_planar_area_m2'],requestedKgPerS=trial,transferredKgPerS=trial,
                            limitedAbsKgPerS=0,pairMassResidualKgPerS=0,volumeSourceM3PerS=trial,
                            phasePowerW=trial*meta['energy_jump_J_kg'],energyJumpJPerKg=meta['energy_jump_J_kg'])
            for phase in names: log+=record('M1C_PHASE_FLUX',time=time,phase=phase,outwardKgPerS=0)
        log+=record('M1_CONTINUUM',solver='compressibleLaserbeamFoam',time=time,dt=meta['delta_t_s'],
                    massKg=2.5e-9,outwardMassFluxKgPerS=0,kineticEnergyJ=0,absorbedLaserPowerW=0,
                    TminK=4101,TmaxK=4101,maxSpeedMS=0.1 if meta['enabled'] else 0)
        direction=-1 if rate<0 else 1
        log+=record('M1C_FLOW',time=time,pressureMeanPa=meta['initial_pressure_Pa']+direction*index*10,
                    pressureMinPa=meta['initial_pressure_Pa'],pressureMaxPa=meta['initial_pressure_Pa']+20)
        log+=record('M1_INTERFACE',time=time,alphaSumMin=1,alphaSumMax=1,analyticMeanAbsError=-1)
        for phase in names:
            mass=2.5e-9-7e-14-transfer if phase=='metal1' else 7e-14+transfer if phase=='metal1vapour' else 0
            log+=record('M1_PHASE',time=time,phase=phase,massKg=mass,alphaVolumeM3=5e-13,nativeSensibleEnergyJ=0)
            log+=record('M1_SPATIAL',time=time,phase=phase,alphaMin=0,alphaMax=1,alphaVolumeM3=5e-13,alphaMomentXM4=0)
    log+='Solving for T, Initial residual = 0, Final residual = 0, No Iterations 0\nEnd\n'
    (case/'log.compressibleLaserbeamFoam').write_text(log)
    (case/'case-result.txt').write_text('exit_status=0\n')
    return log


class InterfaceContract(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.case=Path(self.temp.name)/'case'

    def tearDown(self): self.temp.cleanup()

    def test_integrates_only_final_outer_source(self):
        feedback(self.case)
        result=summary.assess(self.case)
        self.assertTrue(result['passed'],result['gates'])
        self.assertAlmostEqual(result['observations']['transferred_mass_kg'],1.999996e-17,places=25)

    def test_missing_phase_boundary_is_rejected(self):
        log=feedback(self.case)
        log='\n'.join(line for line in log.splitlines() if not (line.startswith('M1C_PHASE_FLUX ') and 'phase=metal1 ' in line))
        (self.case/'log.compressibleLaserbeamFoam').write_text(log)
        with self.assertRaises(ValueError): summary.assess(self.case)

    def test_paired_source_alone_does_not_hide_inventory_error(self):
        log=feedback(self.case)
        # Modify the final phase row without relying on float representation.
        lines=log.splitlines()
        index=max(i for i,line in enumerate(lines) if line.startswith('M1_PHASE ') and 'phase=metal1vapour ' in line)
        row=dict(token.split('=',1) for token in lines[index].split()[1:])
        row['massKg']=float(row['massKg'])+1e-15
        lines[index]=record('M1_PHASE',**row).strip()
        (self.case/'log.compressibleLaserbeamFoam').write_text('\n'.join(lines)+'\n')
        result=summary.assess(self.case)
        self.assertTrue(result['gates']['paired_mass_source'])
        self.assertFalse(result['gates']['vapour_inventory_balance'])

    def test_condensation_has_negative_transfer_and_pressure_response(self):
        feedback(self.case,'condensation')
        result=summary.assess(self.case)
        self.assertTrue(result['passed'],result['gates'])
        self.assertLess(result['observations']['transferred_mass_kg'],0)

    def test_disabled_source_and_generated_open_boundary(self):
        feedback(self.case,'off')
        self.assertTrue(summary.assess(self.case)['passed'])
        open_case=Path(self.temp.name)/'open'
        meta=prepare.prepare(open_case,'open')
        self.assertTrue(meta['open_boundary'])
        self.assertIn('pressureInletOutletVelocity',(open_case/'0/U').read_text())
        self.assertIn('inletValue uniform 1;',(open_case/'0/alpha.metal1vapour').read_text())
        self.assertNotIn('nAlphaSubCycles 2;', (open_case/'system/fvSolution').read_text())


if __name__=='__main__': unittest.main()
