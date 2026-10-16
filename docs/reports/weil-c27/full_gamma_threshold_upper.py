"""Test the cached upper threshold before requesting a common inverse.

Use one existing Arb determinant and the accepted seventeen-positive
inertia bound. A positive upper matrix does not decide the actual
threshold. No Gram, inverse, source action or scientific supplier runs.
"""
import argparse
import hashlib
import platform
from pathlib import Path
import time

import flint
from flint import ctx

from arb_inverse_residual_obstruction import HERE
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, write_packet
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    saved, saved_hash = read_json(args.gram)
    if saved_hash != '0d28fdab5adb68146ae44735809cbfeb2cd0de55bac91143bd4548d5d69501c7':
        raise ValueError('Require the accepted exact-lift finite-Gram packet.')
    names = dict(old='structured-joint-2048-target.json.zlib',
                 generators='joint-residual-generators-2048.json.zlib',
                 pilot='inverse-residual-pilot.json.zlib',
                 base='joint-gamma-6765-rho32-v2-target.json.zlib',
                 boundary='joint-gamma-6765-finite.json.zlib',
                 actions='joint-gamma-6765-rho32-v2-actions.json.zlib',
                 core='joint-gamma-6765-rho32-v2-schur.json',
                 projection='full-gamma-lift-projection.json.zlib',
                 gap='full-gamma-baseline-gap.json.zlib')
    if set(names) != set(saved['packet_sha256']):
        raise ValueError('Incomplete accepted dependency map.')
    guards = {args.gram: saved_hash}
    for key, name in names.items():
        guards[HERE/name] = saved['packet_sha256'][key]
    for key in ('input_sha256', 'helper_sha256'):
        for name, digest in saved[key].items():
            path = HERE/name
            if path in guards and guards[path] != digest:
                raise ValueError('Inconsistent inherited dependency: '+name)
            guards[path] = digest
    source = Path(__file__).resolve()
    for name in ('gamma_boundary_update.py', 'joint_gamma_boundary.py',
                 'structured_inverse_residual_interval.py',
                 'arb_inverse_residual_obstruction.py'):
        path = HERE/name
        guards.setdefault(path, hashlib.sha256(path.read_bytes()).hexdigest())
    guards[source] = hashlib.sha256(source.read_bytes()).hexdigest()
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input/helper hash mismatch: '+path.name)
    old, _ = read_json(HERE/names['old'])
    base, _ = read_json(HERE/names['base'])
    boundary, _ = read_json(HERE/names['boundary'])
    core, _ = read_json(HERE/names['core'])
    projection, _ = read_json(HERE/names['projection'])
    gap, _ = read_json(HERE/names['gap'])
    modes = saved['remaining_modes']
    if (saved['schema'], saved['verdict'], saved['N'], saved['M'],
            saved['finite_first'], saved['finite_last'], saved['grid_bits'],
            len(modes), modes[0], core['remaining_modes'],
            projection['remaining_modes'], gap['remaining_modes'],
            boundary['basis'], core['basis'], projection['basis'],
            len(boundary['bins']), base['metadata']['common_young_parameter'],
            gap['verdict']) != (
            'ACTUAL_C27_LIFT_FINITE_GRAM_V1', 'same-lift-common-finite-gram-enclosed',
            1280, 32768, 2049, 32768, 640, 18, 0, modes, modes, modes,
            old['metadata']['basis'], old['metadata']['basis'], old['metadata']['basis'],
            19, '1/32', 'same-Wbar-stronger-base-fixed-dual-update-excluded'):
        raise ValueError('Different support, basis, lift order, channels or far applicability.')
    if (projection['original_metadata_grid_bits'],
            base['metadata']['original_metadata_grid_bits']) != (512, 512):
        raise ValueError('Keep allowance metadata and entry export grids distinct.')
    if saved['source_sha256'] != saved['helper_sha256']['full_gamma_lift_finite_gram.py']:
        raise ValueError('Different finite-Gram producer.')
    rows = saved['matrices']['threshold_upper']
    if len(rows) != 18 or any(len(row) != 18 for row in rows):
        raise ValueError('Require all eighteen common coordinates.')
    if any(int(rows[i][j][0]) > int(rows[i][j][1]) or rows[i][j] != rows[j][i]
           for i in range(18) for j in range(18)):
        raise ValueError('Require valid symmetric entry enclosures, without midpoint replacement.')
    upper = matrix_in(rows, saved['grid_bits'])
    if not all(value.is_finite() for value in upper.entries()):
        raise ValueError('Nonfinite cached upper matrix.')
    determinant = upper.det()
    # Reuse the degree-aware absolute grid of small_schur_rank_floor.py.
    # A degree-eighteen determinant must not be rounded on the entry grid.
    determinant_bits = 18*saved['grid_bits']+64
    pair = envelope(determinant, determinant_bits)
    sign = (1 if pair is not None and int(pair[0]) > 0 else
            -1 if pair is not None and int(pair[1]) < 0 else 0)
    verdict = {1: 'upper-positive-actual-threshold-undecided',
               -1: 'upper-negative-direction-no-strict-negative-lower-source',
               0: 'inconclusive-upper-determinant'}[sign]
    result = dict(schema='ACTUAL_C27_THRESHOLD_UPPER_DETERMINANT_V1',
                  finite_gram_sha256=saved_hash,
                  inherited_packet_sha256=saved['packet_sha256'],
                  input_sha256=saved['input_sha256'], remaining_modes=modes,
                  source_sha256=guards[source],
                  consumed_sha256={path.name: digest for path, digest in guards.items()},
                  working_bits=args.precision, matrix_entry_grid_bits=saved['grid_bits'],
                  determinant_grid_bits=determinant_bits, determinant_interval=pair,
                  determinant_sign=sign, native_determinant_calls=1,
                  runtime=dict(python=platform.python_version(), flint=flint.__version__,
                               machine=platform.machine()),
                  verdict=verdict, elapsed_seconds=time.monotonic()-started,
                  conditions=list(dict.fromkeys(saved['conditions']+gap['conditions']))+[
                      'The accepted same-object channel/inverse comparison gives N>=0 and 0<=Theta<=SJ.',
                      'Accepted Xi has exactly17positive1negative; SJ>=0 implies U=Xi+SJ has at least17positives.',
                      'The cached entry enclosure contains the same exact U, not an entry-endpoint Loewner bound.',
                      'The accepted finite source identity n_-(Csharp)=n_+(Xi+Theta)-17 applies unchanged.'],
                  scope='One new cached18x18 upper determinant only. A positive upper does not settle '
                    'Xi+Theta or certify a negative source. A negative upper direction only excludes '
                    'strict negative directions of the same finite-J/fixed18lift lower source comparison; '
                    'it does not exclude actual sufficient updates generally or settle endpoint cases. '
                    'A zero-containing/nonfinite determinant is inconclusive for this evaluation. '
                    'No old Gram, source action, inverse/solve, scientific supplier, generic theorem, '
                    'algorithm, Lean, actual retained/cofinal positivity, Robin-tail or RH result.')
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('A consumed input/helper changed: '+path.name)
    write_packet(args.out, result)
    print('determinant', determinant.str(18), flush=True)
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gram', type=Path, default=HERE/'full-gamma-lift-finite-gram.json.zlib')
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib' or args.precision < 832:
        parser.error('Use a fresh .zlib output and at least192 guardbits over the inherited entry grid.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
