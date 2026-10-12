"""Reference consistency and diagnostic rejection; not native CFD evidence."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPTS=Path(__file__).parents[1]/'scripts'
def load(name,file):
    spec=importlib.util.spec_from_file_location(name,SCRIPTS/file)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
reference=load('spatial_oracle','audit-m2a-spatial-consistency.py')
prepare=load('spatial_prepare','prepare-m2a-joint.py')
summary=load('spatial_summary','summarize-m2a-spatial.py')


class SpatialTests(unittest.TestCase):
    def test_continuum_sampling_has_nonzero_wall_defect_and_refines(self):
        for mode in ('fixed','rotate'):
            cases=[reference.fixture(mode,n,0) for n in (36,48,64)]
            self.assertTrue(all(c['wall_rms_m_s']>1e-7 for c in cases))
            self.assertGreater(cases[0]['wall_rms_m_s'],cases[-1]['wall_rms_m_s'])
            self.assertGreater(abs(cases[0]['stress_load_ratio']-1),abs(cases[-1]['stress_load_ratio']-1))
            self.assertLess(abs(cases[-1]['stress_load_ratio']-1),.06)

    def test_no_advance_reference_and_missing_or_corrupt_native_evidence(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            for name in prepare.CASES:
                if name in ('fixed-restart','fixed-mpi2'):continue
                case=root/name;ranks=8 if 'finer' in name else 1
                meta=prepare.prepare(case,name,ranks,spatial_only=True);r=meta['spatial_reference']
                self.assertIn('endTime 0;',(case/'system/controlDict').read_text())
                self.assertIn('spatialAuditOnly true;',(case/'constant/mechanicalProperties').read_text())
                (case/'result.txt').write_text('exit_status=0\n')
                (case/'log.solver').write_text(f'nProcs : {ranks}\nM2A_SPATIAL_REFERENCE markers={r["markers"]} wallRms={r["wall_rms_m_s"]} wallMax={r["wall_max_m_s"]} stressRatio={r["stress_load_ratio"]}\nM2A_EXECUTION_COMPLETE\n')
            self.assertTrue(summary.summarize(root)['passed'])
            case=root/'fixed-finer';log=case/'log.solver';text=log.read_text()
            for bad in (text.replace('nProcs : 8','nProcs : 1'),text.replace('M2A_SPATIAL_REFERENCE','M2A_OLD_REFERENCE'),text.replace('stressRatio=','stressRatio=nan #')):
                log.write_text(bad);self.assertFalse(summary.summarize(root)['passed'])
                json.dumps(summary.summarize(root),allow_nan=False)
            with self.assertRaises(ValueError):prepare.prepare(root/'bad','fixed-restart',1,True)


if __name__=='__main__':unittest.main()
