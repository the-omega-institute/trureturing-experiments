"""Apply classical Sylvester's criterion to a saved small Schur enclosure.

Reuse the certified positive prefix and the validated dense complement.
FLINT's existing ball determinant supplies each leading principal minor;
no target production, solve, eigensolver or displacement elimination runs.
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
    target, prefix, remaining, target_hash, pivot_raw = load_target(
        args.target, args.pivots, args.max_remainder)
    core, core_hash = read_json(args.schur)
    if (core['schema'], core['target_packet_sha256'], core['pivot_result_sha256'],
            core['remaining_modes'], core['reused_positive_coordinate_count'],
            core['input_sha256'], core['N'], core['M'], core['last_mode'], core['basis']) != (
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1', target_hash,
            hashlib.sha256(pivot_raw).hexdigest(), remaining, len(prefix),
            target['metadata']['input_sha256'], target['metadata']['N'],
            target['metadata']['M'], target['metadata']['last_mode'], target['metadata']['basis']):
        raise ValueError('Require the same complete target, positive prefix and dense Schur complement.')
    matrix = matrix_in(core['reduced_schur_intervals'], core['grid_bits'])
    size = len(remaining)
    if (matrix.nrows(), matrix.ncols()) != (size, size):
        raise ValueError('Incomplete saved small complement.')
    minors, verdict = [], 'positive-leading-minors'
    for k in range(1, size+1):
        minor = arb_mat([[matrix[i, j] for j in range(k)] for i in range(k)]).det()
        # The minors can be far smaller than individual entries. Keep a
        # separate absolute grid for each degree instead of rounding to0.
        bits = k*core['grid_bits']+64
        pair = envelope(minor, bits)
        sign = ('positive' if pair is not None and int(pair[0]) > 0 else
                'negative' if pair is not None and int(pair[1]) < 0 else 'inconclusive')
        minors.append(dict(order=k, grid_bits=bits, determinant_interval=pair, sign=sign))
        print('leading minor', k, sign, minor.str(16), flush=True)
        if sign != 'positive':
            verdict = 'negative-leading-minor' if sign == 'negative' else 'inconclusive-leading-minor'
            break
    result = dict(schema='ACTUAL_C27_SMALL_SCHUR_SYLVESTER_V1',
                  target_packet_sha256=target_hash,
                  pivot_result_sha256=hashlib.sha256(pivot_raw).hexdigest(),
                  schur_result_sha256=core_hash,
                  reused_positive_coordinate_count=len(prefix), remaining_modes=remaining,
                  working_bits=args.precision, verdict=verdict, leading_minors=minors,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  conditions=list(dict.fromkeys(core['conditions']+json.loads(pivot_raw).get('conditions', [])))+[
                      'The saved real symmetric Schur enclosure contains the exact complement of this target.',
                      'All strict leading-minor signs are read from outward endpoint pairs.'],
                  scope='Classical Sylvester sign test of one saved small complement. '
                    'All positive leading minors certify this sufficient finite target under inherited premises. '
                    'A negative leading minor rejects only that lower target; inconclusive decides neither sign. '
                    'No new generic method, Lean verification, cofinal positivity or RH result.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(verdict, 'seconds', time.monotonic()-started, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('target', 'pivots', 'schur', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--precision', type=int, default=1536)
    parser.add_argument('--max-remainder', type=int, default=32)
    args = parser.parse_args()
    if args.out.exists() or args.precision < 832:
        parser.error('Use a fresh output and at least832 working bits.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
