"""Reproduce X sections 9--10 without rewriting the published coefficient data."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

import flint
import sympy


ROOT = Path(__file__).resolve().parent
HASHES = {
    'R.json': 'a7b5c02552dea6e43bbd4919bb0f0db88e92e838fca26b1bbc976603437d568e',
    'E.json': 'de87cab361a10a0f7a84a15a3738b5e8cad7702331698b6fcf2f3b1e530aff36',
    'H_factor_1.json': 'a4c611369cd5eedc73c7314581a507da880b7db30f5d2fc6fd01eb9866d5229e',
    'collision_support.json': 'ce0e91ac8fcf7b8637f0520f368938a9f165d305cc76c156c949cce5bbb727f1',
    'support_isolations.json': '9f0cd27bde73dcbedac70eca973ee9942ee7d8bc46998232ec168e471ab62cd7',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def isolation_data(rows):
    return [{k: v for k, v in row.items() if k != 'seconds'} for row in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('document', type=Path, help='the X theory volume containing sections 8--10')
    parser.add_argument('--output', type=Path, help='optional compact mathematical result JSON')
    args = parser.parse_args()
    assert sys.version_info >= (3, 12)
    document = args.document.resolve()
    advertised = ROOT / 'global_preflight'
    assert all(digest(advertised / name) == value for name, value in HASHES.items())
    with TemporaryDirectory(prefix='.reproduce-', dir=ROOT) as temporary:
        work = Path(temporary)
        for name in ('global_preflight', 'bprime_display_correction', 'dependencies'):
            shutil.copytree(ROOT / name, work / name)
        shutil.copyfile(ROOT / 'scan_fibers.py', work / 'scan_fibers.py')
        commands = [
            ['global_preflight/probe.py'],
            ['global_preflight/discriminant.py'],
            ['global_preflight/support.py'],
            ['global_preflight/isolate_support.py'],
            ['global_preflight/verify_certificate.py'],
            ['bprime_display_correction/verify_bprime.py', str(document), 'bprime_results.json'],
        ]
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        for command in commands:
            completed = subprocess.run([sys.executable, *command], cwd=work, env=env,
                                       text=True, capture_output=True)
            print(json.dumps({'program': command[0], 'exit': completed.returncode}), flush=True)
            if completed.returncode:
                raise RuntimeError(completed.stdout + completed.stderr)
        generated = work / 'global_preflight'
        regenerated = [name for name in HASHES if name != 'support_isolations.json']
        assert all(digest(generated / name) == HASHES[name] for name in regenerated)
        isolated = read(generated / 'support_isolations.json')
        assert isolation_data(isolated) == isolation_data(read(advertised / 'support_isolations.json'))
        counts = [row['count'] for row in read(generated / 'support_results.json')
                  if row['label'] == 'factor_real_root_count']
        assert counts == [row['real_roots'] for row in isolated] == [1, 2, 2, 2, 4, 20, 24, 8, 49]
        assert [row['degree'] for row in isolated] == [1, 2, 2, 2, 4, 24, 32, 44, 91]
        identities = {k: v for k, v in read(generated / 'verification_results.json').items()
                      if k != 'elapsed_seconds'}
        bprime = read(work / 'bprime_results.json')
        assert all(identities.values()) and bprime['passed']
        report = {
            'passed': True,
            'requirements': {'python': sys.version.split()[0], 'sympy': sympy.__version__,
                             'python-flint': flint.__version__},
            'advertised_sha256': HASHES,
            'regenerated_coefficient_files': regenerated,
            'regenerated_isolations_match_except_timing': True,
            'global_exact_identities': identities,
            'support_degree': 202,
            'support_factor_degrees': [row['degree'] for row in isolated],
            'distinct_real_roots_per_factor': counts,
            'distinct_real_candidate_responses': sum(counts),
            'bprime_checks': len(bprime['checks']),
            'bprime_checks_passed': all(row['passed'] for row in bprime['checks']),
            'bprime_measurements': bprime['measurements'],
            'arithmetic': 'rational polynomial identities, exact Sturm counts and rational root isolation; outward integer Bprime intervals at scale 10^80',
            'scope': 'fixed parent (3/7,-35/29); X section 9 full Bprime band and section 10 necessary MM support',
            'not_established': ['112 physical specializations', 'C_MM intersection C_MJ emptiness',
                                'root eight/nine-read decision', 'nine-read upper bound',
                                'Lean verification', 'independent review'],
        }
    assert all(digest(advertised / name) == value for name, value in HASHES.items())
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'global_identity_checks': len(identities),
                      'bprime_checks': report['bprime_checks'], 'necessary_candidates': sum(counts)}))


if __name__ == '__main__':
    main()
