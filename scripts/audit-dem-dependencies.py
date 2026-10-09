"""User-side read-only inventory. Does not execute or build OpenFOAM/DEM."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
from datetime import datetime


def inventory(roots, max_depth=5):
    headers, libraries, executables, errors = [], [], [], []
    skip = {'.git', '.cache', '.local', 'node_modules', '.runs', '.build', 'proc', 'sys', 'dev'}
    for root in roots:
        root = Path(root).expanduser().absolute()
        if not root.is_dir():
            continue
        def onerror(e):
            errors.append(str(e))
        for folder, dirs, files in os.walk(root, followlinks=False, onerror=onerror):
            depth = len(Path(folder).relative_to(root).parts)
            dirs[:] = [] if depth >= max_depth else [d for d in dirs if d not in skip and not d.startswith('.')]
            for name in files:
                p = Path(folder) / name
                if name == 'library.h':
                    try:
                        data = p.read_bytes()
                        if b'lammps_open' in data and b'lammps_gather_atoms' in data:
                            headers.append({'path': str(p), 'source_directory': str(p.parent),
                                            'sha256': hashlib.sha256(data).hexdigest()})
                    except OSError as e:
                        errors.append(str(e))
                if (name.startswith(('libliggghts', 'liblammps', 'liblmp')) and '.so' in name):
                    libraries.append({'path': str(p), 'resolved_path': str(p.resolve()),
                                      'exists': p.exists()})
                if name in {'liggghts', 'lmp_auto', 'lmp_mpi'}:
                    executables.append(str(p))
    return {'searched_roots': [str(p) for p in roots], 'max_depth': max_depth,
            'headers': headers, 'shared_libraries': libraries,
            'executable_candidates': executables, 'read_errors': errors,
            'scope': 'Inventory only; does not establish header/library match, MPI ABI or execution feasibility.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', action='append', type=Path,
                        help='Additional installation root; default search covers HOME and /opt.')
    parser.add_argument('--depth', type=int, default=5)
    args = parser.parse_args()
    if not 1 <= args.depth <= 10:
        parser.error('depth must be between 1 and 10')
    repo = Path(__file__).resolve().parents[1]
    report = repo / '.runs' / ('dem-dependency-audit-' + datetime.now().strftime('%Y%m%d-%H%M%S') + '-' + str(os.getpid()))
    report.mkdir(parents=True)
    roots = [Path.home(), Path('/opt')] + (args.root or [])
    data = inventory(roots, args.depth)
    data['explicit_environment'] = {k: os.environ.get(k) for k in
                                   ('LIGGGHTS_SRC', 'LIGGGHTS_LIB', 'WM_PROJECT_VERSION', 'WM_OPTIONS', 'WM_MPLIB')}
    commit = subprocess.run(['git', '-C', str(repo), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    (report / 'commit.txt').write_text(commit.stdout or commit.stderr)
    (report / 'dependency-inventory.json').write_text(json.dumps(data, indent=2) + '\n')
    # No shell/alias invocation, binary execution, recursive symlink walk or installation changes.
    (report / 'result.txt').write_text('last_stage=inventory\nexit_status=0\nvalidation=not_run\n')
    archive = report.with_name(report.name + '.tar.gz')
    with tarfile.open(archive, 'w:gz') as t:
        t.add(report, arcname=report.name)
    print('Candidate headers:')
    for item in data['headers']:
        print('  ' + item['path'])
    print('Candidate shared libraries:')
    for item in data['shared_libraries']:
        print('  ' + item['path'])
    if not data['headers'] or not data['shared_libraries']:
        print('Required files not both found within search scope; absence here is not proof of no installation.')
    print('Feedback archive: ' + str(archive))


if __name__ == '__main__':
    main()
