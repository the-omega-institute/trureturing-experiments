"""Portable consumer of the executed fixed-source interval routines."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--verify-files', action='store_true')
    mode.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('Use ordinary Python without -O.')
    manifest = json.loads((ROOT / 'artifact-manifest.json').read_text())
    for entry in manifest['files']:
        raw = (ROOT / entry['path']).read_bytes()
        if len(raw) != entry['bytes'] or digest(raw) != entry['sha256']:
            raise RuntimeError('Package binding mismatch: ' + entry['path'])
    if args.verify_files:
        print(json.dumps(dict(files_verified=len(manifest['files']),
                              science_arithmetic_run=False)))
        return

    def load(relative, name):
        path = ROOT / relative
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module

    core = load('r2380/diagonal_bound_query.py', 'portable_core')
    base = load('r2429/signed_interval_replay.py', 'portable_base')
    small = load('r2490/small_boundary_inertia.py', 'portable_small')
    fixed = load('r2510/fixed_congruence.py', 'portable_fixed')
    contraction_raw = (ROOT / 'signed-eta1-contractions-result-r2473.json').read_bytes()
    reduction_raw = (ROOT / 'signed-eta1-border-reduction-result-r2499.json').read_bytes()
    packet = json.loads(reduction_raw)
    if digest(reduction_raw) != fixed.REDUCTION_SHA256:
        raise RuntimeError('Original reduction mismatch.')
    recomputed = small.reduce_paid_boundary(core, base, contraction_raw,
                                            digest(contraction_raw))
    for key, value in recomputed.items():
        if packet.get(key) != value:
            raise RuntimeError('Recomputed reduction differs at ' + key)
    basis_raw = (ROOT / 'actual-fixed-dyadic-basis.json').read_bytes()
    basis = json.loads(basis_raw)
    saved = json.loads((ROOT / 'signed-eta1-preconditioned112-result-r2512.json').read_text())
    if digest(basis_raw) != saved['certificate_bindings']['basis_sha256']:
        raise RuntimeError('Executed basis mismatch.')
    order = basis['permutation']
    blocks = basis['original_label_blocks']
    if (basis['dimension'] != 112 or basis['basis_grid_bits'] != 4096
            or basis['original_arithmetic_bits'] != 1024
            or sorted(order) != list(range(112))
            or [i for block in blocks for i in block] != order):
        raise RuntimeError('Invalid exact basis metadata.')
    small._check_blocks(core, order, blocks, complete=True)
    original = small._decode_upper(core, packet['reduced112_upper'], 112, 1024,
                                   'original112')
    permuted = [[original[i][j] for j in order] for i in order]
    matrix = fixed.enclose_congruence(core, base, permuted,
                                      basis['basis_upper_point_integers'], 4096)
    _, active, inertia, transcript = small._eliminate(core, base, matrix, order, blocks)
    core.require(not active and sum(inertia) == 112, 'incomplete_coverage',
                 'No count from an incomplete suffix.')
    core.require(inertia[1] >= 56, 'inherited_contract_inconsistency',
                 'Restricted negative index cannot be negative.')
    if inertia != saved['original_reduced_inertia'] or transcript != saved['pivot_transcript']:
        raise RuntimeError('Complete replay differs from saved execution.')
    print(json.dumps(dict(original_reduced_inertia=inertia,
                          border_inertia=[inertia[0], inertia[1]+21, 0],
                          same_source_restricted_inertia=[1205+inertia[0], inertia[1]-56, 0],
                          inherited_source_and_action_contracts_required=True,
                          eta_budget_or_RH_claim=False)))


if __name__ == '__main__':
    main()
