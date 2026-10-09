"""Offline tests of feedback handling, not CFD validation."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('portability',ROOT/'scripts/summarize-m1b-portability.py')
p=importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


class PortabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)

    def tearDown(self): self.temp.cleanup()

    def case(self,name,temperature=4101,energy_increment=1e-9):
        case=self.root/name; case.mkdir()
        (case/'m1b-case.json').write_text(json.dumps(dict(length_m=1e-4,vapour_reference_energy_offset_J_per_kg=0)))
        lines=[]
        for i in range(3):
            lines.append(f'M1_CONTINUUM time={i*1e-9} massKg=2e-9 TminK=4101 TmaxK={temperature} kineticEnergyJ=0 absorbedLaserPowerW=1')
            lines.append(f'M1_PHASE time={i*1e-9} phase=metal1vapour nativeSensibleEnergyJ={0.01+i*energy_increment} massKg=1e-12')
            lines.append(f'M1_SPATIAL time={i*1e-9} phase=metal1vapour alphaVolumeM3=5e-13')
        (case/'log.compressibleLaserbeamFoam').write_text('\n'.join(lines))
        return case

    def test_equal_histories(self):
        self.assertTrue(p.compare(self.case('base'),self.case('other'))['passed'])

    def test_temperature_mismatch_rejected(self):
        self.assertFalse(p.compare(self.case('base'),self.case('other',temperature=4101.001))['passed'])

    def test_energy_change_uses_signal_not_large_initial_energy(self):
        result=p.compare(self.case('base'),self.case('other',energy_increment=1.001e-9))
        self.assertFalse(result['gates']['energy_change_absolute_J'])

    def test_merge_retains_steps_and_rejects_bad_checkpoint(self):
        first=self.root/'first'; second=self.root/'second'
        first.write_text('M1_CONTINUUM time=0\nM1_CONTINUUM time=1\nM1_PHASE time=1 phase=x\nEnd\n')
        second.write_text('M1_CONTINUUM time=1\nM1_PHASE time=1 phase=x\nM1_CONTINUUM time=2\nM1_PHASE time=2 phase=x\nEnd\n')
        merged=p.merge(first,second)
        self.assertEqual([r['time'] for r in p.integrated.rows(merged,'M1_CONTINUUM')],['0','1','2'])
        self.assertEqual(merged.splitlines().count('End'),1)
        second.write_text('M1_CONTINUUM time=0\nM1_CONTINUUM time=2\nEnd\n')
        with self.assertRaises(ValueError): p.merge(first,second)
        second.write_text('M1_CONTINUUM time=1\nM1_CONTINUUM time=2\n')
        with self.assertRaises(ValueError): p.merge(first,second)


if __name__=='__main__': unittest.main()
