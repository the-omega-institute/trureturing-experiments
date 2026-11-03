"""Immutable public source and existing exact experiment suppliers."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import subprocess
import sys
import types

ROOT = Path(__file__).resolve().parent.parent
THEORY_REPO = Path('/Users/auricstudio/trureturing')
THEORY_PIN = '0434b0714812f29d6799f8e4ebd958f5eab2f987'
BASE = 'c27e6de785aefbf37c0789a3f4cf63c6667ff3d7'
EXPERIMENT_PIN = 'bfd14367a5cbab3ee85cc0ce13db70dd797213b3'
SUPPLIER_PIN = '865f38257ead1b4e857ff915c23b0b04276bd7f3'
POSITIVE_PIN = 'b490f4f17c97645ce7c6495ca8635bce3d3a1af7'
THEORY_PATH = 'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_MONOTONE_LAB_INDEXED_EXECUTION.md'
THEORY_HASH = '0f93b9f2634965902d7ea35ee675e1886d2309ec10b017e926ff756f3136fcd2'
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
    module.__file__ = f'{pin}:{path}'
    sys.modules[name] = module
    exec(compile(source, module.__file__, 'exec'), module.__dict__)
    return module


def suppliers():
    # Exact older inversion supplier is used only at <=4 *active* indices.
    old = load(SUPPLIER_PIN, 'rro-structural-recovery/exact_checks.py', 'rro_old')
    # The general evaluator imports this exact dependency by its original name.
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


def source_manifest():
    blob(THEORY_REPO, THEORY_PIN, THEORY_PATH, THEORY_HASH)
    blob(THEORY_REPO, THEORY_PIN, 'Library/Computability/nearywoodsmurphyglaschick2014wang.md',
         '53d0eb487e0de645e05da928fdae12a2c4abc70b8f784e3f062afead786357a1')
    for name in ('ATOMIC_GENERATION_ACQUISITION', 'UNIFIED_JOINT_LAW_BRIDGE', 'PROCESS_GEOMETRY'):
        blob(THEORY_REPO, BASE, f'docs/develop/theory/RECURSIVE_RELATIONAL_OBSERVATION_{name}.md')
    return SOURCES
