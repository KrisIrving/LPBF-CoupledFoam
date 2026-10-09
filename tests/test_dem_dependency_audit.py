import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('audit', Path(__file__).resolve().parents[1]
                                            / 'scripts/audit-dem-dependencies.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class DependencyInventory(unittest.TestCase):
    def test_detects_interface_not_arbitrary_library_header_and_honours_depth(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'src').mkdir()
            (root / 'src/library.h').write_text('lammps_open lammps_gather_atoms')
            (root / 'lib').mkdir()
            (root / 'lib/libliggghts.so').write_text('fixture only')
            (root / 'other').mkdir()
            (root / 'other/library.h').write_text('unrelated library')
            deep = root / 'a/b/c'
            deep.mkdir(parents=True)
            (deep / 'library.h').write_text('lammps_open lammps_gather_atoms')
            result = audit.inventory([root], max_depth=2)
            self.assertEqual(len(result['headers']), 1)
            self.assertEqual(len(result['shared_libraries']), 1)
            self.assertEqual((root / 'lib/libliggghts.so').read_text(), 'fixture only')

    def test_missing_root_and_excluded_reports_do_not_create_candidates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / '.runs').mkdir()
            (root / '.runs/library.h').write_text('lammps_open lammps_gather_atoms')
            result = audit.inventory([root, root / 'missing'])
            self.assertFalse(result['headers'])
            self.assertFalse(result['shared_libraries'])


if __name__ == '__main__':
    unittest.main()
