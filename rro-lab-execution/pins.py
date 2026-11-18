"""Candidate content binding and immutable mathematical/experiment suppliers."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import types

ROOT = Path(__file__).resolve().parent.parent
THEORY_PIN = '0434b0714812f29d6799f8e4ebd958f5eab2f987'
BASE = 'c27e6de785aefbf37c0789a3f4cf63c6667ff3d7'
EXPERIMENT_PIN = 'bfd14367a5cbab3ee85cc0ce13db70dd797213b3'
SUPPLIER_PIN = '865f38257ead1b4e857ff915c23b0b04276bd7f3'
POSITIVE_PIN = 'b490f4f17c97645ce7c6495ca8635bce3d3a1af7'
THEORY_PATH = 'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_MONOTONE_LAB_INDEXED_EXECUTION.md'
BIBLIOGRAPHY_PATH = 'Library/Computability/nearywoodsmurphyglaschick2014wang.md'
SOURCES = {}


def blob(repo, pin, path, expected=None):
    descriptors = tuple(int(x) for x in os.environ.get('WORKTREE_SCOPE_FDS','').split(',') if x.isdigit())
    data = subprocess.check_output(['git', '-C', str(repo), 'show', f'{pin}:{path}'],pass_fds=descriptors)
    digest = hashlib.sha256(data).hexdigest()
    if expected is not None:
        assert digest == expected, (path, digest)
    SOURCES[f'{pin}:{path}'] = {'sha256': digest, 'bytes': len(data)}
    return data


def load(pin, path, name):
    source = blob(ROOT, pin, path)
    module = types.ModuleType(name)
    # Legacy suppliers derive their Git root from __file__; keep it absolute.
    # The compiled source itself remains the immutable Git blob.
    module.__file__ = str(ROOT/path)
    sys.modules[name] = module
    exec(compile(source, f'{pin}:{path}', 'exec'), module.__dict__)
    return module


def suppliers():
    old = load(SUPPLIER_PIN, 'rro-structural-recovery/exact_checks.py', 'rro_old')
    saved = sys.modules.get('exact_checks')
    sys.modules['exact_checks'] = old
    try:
        evaluator = load(EXPERIMENT_PIN, 'rro-structural-recovery/exact_evaluator.py', 'rro_general')
        positive = load(POSITIVE_PIN, 'rro-positive-certificates/exact_checks.py', 'rro_positive')
    finally:
        if saved is None:
            sys.modules.pop('exact_checks', None)
        else:
            sys.modules['exact_checks'] = saved
    return old, evaluator, positive


def source_manifest(theory_root, candidate_manifest):
    """Fail closed on changed candidate bytes; no fallback to a sealed old draft."""
    SOURCES.clear()
    theory_root, candidate_manifest = Path(theory_root), Path(candidate_manifest)
    manifest_bytes = candidate_manifest.read_bytes()
    candidate = json.loads(manifest_bytes)
    assert candidate['schema'] == 1 and candidate['sealed_unit1'] == THEORY_PIN
    assert set(candidate['files']) == {THEORY_PATH, BIBLIOGRAPHY_PATH}
    for name, expected in candidate['files'].items():
        data = (theory_root/name).read_bytes()
        observed = {'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}
        assert observed == expected, ('Candidate content mismatch', name, observed)
        SOURCES['candidate:'+name] = observed
    SOURCES['candidate-manifest'] = {'sha256': hashlib.sha256(manifest_bytes).hexdigest(),
                                     'bytes': len(manifest_bytes)}
    for name in ('ATOMIC_GENERATION_ACQUISITION', 'UNIFIED_JOINT_LAW_BRIDGE', 'PROCESS_GEOMETRY'):
        blob(theory_root, BASE, f'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_{name}.md')
    return SOURCES
