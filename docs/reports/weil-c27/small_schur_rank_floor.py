"""Apply Jacobi's signature criterion to a saved small Schur complement.

FLINT evaluates leading determinants; their sign changes count negative
pivots. Stop once a requested rank of positive updates is insufficient,
or a determinant is inconclusive. No matrix construction or solve runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import arb_mat, ctx

from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope
from structured_residual_schur_finish import load_target


def consume(args):
    started = time.monotonic()
    target, prefix, remaining, target_hash, pivot_raw = load_target(args.target, args.pivots)
    core, core_hash = read_json(args.schur)
    if (core['schema'], core['target_packet_sha256'], core['pivot_result_sha256'],
            core['remaining_modes'], core['reused_positive_coordinate_count'], core['input_sha256']) != (
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1', target_hash,
            hashlib.sha256(pivot_raw).hexdigest(), remaining, len(prefix),
            target['metadata']['input_sha256']):
        raise ValueError('Require the same complete base, positive prefix and saved complement.')
    matrix = matrix_in(core['reduced_schur_intervals'], core['grid_bits'])
    if (matrix.nrows(), matrix.ncols()) != (len(remaining), len(remaining)):
        raise ValueError('Incomplete small complement.')
    previous, negative = 1, 0
    minors, verdict = [], 'full-complement-inertia-counted'
    for k in range(1, len(remaining)+1):
        determinant = arb_mat([[matrix[i, j] for j in range(k)] for i in range(k)]).det()
        bits = k*core['grid_bits']+64
        pair = envelope(determinant, bits)
        sign = 1 if pair is not None and int(pair[0]) > 0 else (
            -1 if pair is not None and int(pair[1]) < 0 else 0)
        pivot_sign = sign*previous
        if sign:
            negative += int(pivot_sign < 0)
        minors.append(dict(order=k, grid_bits=bits, determinant_interval=pair,
                           determinant_sign=sign, pivot_sign=pivot_sign,
                           negative_count_through_order=negative))
        print('leading order', k, 'det sign', sign, 'negative count', negative, flush=True)
        if not sign:
            verdict = 'inconclusive-leading-determinant'
            break
        previous = sign
        if negative > args.update_rank:
            verdict = 'positive-update-rank-insufficient'
            break
    result = dict(schema='ACTUAL_C27_SMALL_SCHUR_RANK_FLOOR_V1',
                  target_packet_sha256=target_hash,
                  pivot_result_sha256=hashlib.sha256(pivot_raw).hexdigest(),
                  schur_result_sha256=core_hash, remaining_modes=remaining,
                  reused_positive_coordinate_count=len(prefix), working_bits=args.precision,
                  tested_update_rank=args.update_rank, verdict=verdict,
                  negative_count_lower=negative, leading_minors=minors,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  conditions=list(dict.fromkeys(core['conditions']+json.loads(pivot_raw).get('conditions', [])))+[
                      'The saved real symmetric enclosure contains the exact complement of this base.',
                      'Nonzero leading determinant signs determine the principal-block inertia by Jacobi signature.',
                      'A positive-semidefinite update of rank r cannot remove more than r negative directions.'],
                  scope='A lower bound on negative inertia of this sufficient lower base only. '
                    'If the count exceeds r, no positive update of rank at most r makes that base PSD. '
                    'This does not exclude a larger-rank gain or additional comparison slack, and '
                    'does not give a negative direction of the actual Schur form. No full matrix '
                    'solve, spectral retry, new generic method, Lean verification or RH result.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(verdict, 'seconds', time.monotonic()-started, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('target', 'pivots', 'schur', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--precision', type=int, default=1536)
    parser.add_argument('--update-rank', type=int, default=1)
    args = parser.parse_args()
    if args.out.exists() or args.precision < 832 or args.update_rank < 0:
        parser.error('Use a fresh output, at least832 working bits and nonnegative update rank.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
