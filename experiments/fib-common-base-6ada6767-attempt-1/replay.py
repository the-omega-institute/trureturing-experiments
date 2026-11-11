"""Replay retained exact certificates while denying reads from provenance checkouts.

This wrapper checks file retention and layout only; all mathematical checks remain
in their original programs. It runs no numerical discovery or root solver.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import runpy
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--theory', type=Path, required=True)
args = parser.parse_args()
theory = args.theory.resolve(strict=True)
inventory = json.loads((HERE / 'DEPENDENCIES.json').read_text())
legacy_roots = [Path(s['checkout']).resolve() for s in inventory['origins'].values()
                if s['kind'] == 'read-only experiment source']
local_reads = set()


def audit(event, values):
    if event != 'open' or not isinstance(values[0], (str, bytes, os.PathLike)):
        return
    path = Path(os.fsdecode(values[0])).resolve()
    if any(path == old or old in path.parents for old in legacy_roots):
        raise RuntimeError('provenance checkout access forbidden: ' + str(path))
    if ROOT in path.parents:
        local_reads.add(str(path.relative_to(ROOT)))


sys.addaudithook(audit)
if sys.flags.optimize:
    parser.error('assertions are required; do not use -O or PYTHONOPTIMIZE')
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))


def check_inventory():
    for item in inventory['retained_files'] + inventory['added_files']:
        path = ROOT / item['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'], item['path']


check_inventory()
local_reads.clear()
# Run the exact entry points, preserving their CLI and mathematical bodies.
for name, extra in [('reuse_inverse.py', []),
                    ('verify_exact.py', ['--theory', str(theory)]),
                    ('verify_paper_certificate.py', ['--theory', str(theory)])]:
    sys.argv = [str(HERE / name)] + extra
    runpy.run_path(sys.argv[0], run_name='__main__')

# Check the two discovery-only supplier edges without starting a search.
sys.argv = [str(HERE / 'discover.py')]
for consumer, variable in [('discover.py', 'REUSED'), ('joint_search.py', 'code')]:
    tree = ast.parse((HERE / consumer).read_text())
    node = next(n for n in ast.walk(tree) if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == variable for t in n.targets))
    value = eval(compile(ast.Expression(node.value), consumer, 'eval'),
                 {'Path': Path, '__file__': str(HERE / consumer), 'ast': ast})
    if variable == 'REUSED':
        runpy.run_path(str(value), run_name='__main__')
    else:
        blocks = next(ast.literal_eval(n.value) for n in value.body
                      if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == 'blocks' for t in n.targets))
        assert len(blocks) == 3
supplier_reads = sorted(p for p in local_reads
                        if p.startswith('experiments/fib-self-calibration/'))
expected_suppliers = {edge['supplier'] for edge in inventory['supplier_edges']
                      if edge['supplier'].startswith('experiments/fib-self-calibration/')}
assert expected_suppliers == set(supplier_reads)
check_inventory()
print('PORTABLE_REPLAY_PASS; retained hashes unchanged; provenance checkout reads forbidden')
print('LOCAL_SUPPLIER_READS', json.dumps(supplier_reads))
