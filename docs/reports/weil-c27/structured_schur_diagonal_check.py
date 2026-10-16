"""Read saved Schur enclosures and check a negative diagonal using integers.

This consumes the complete target, positive-prefix certificate and validated
dense-Schur result. It does not produce a target, solve or spectrum. The
conclusion is conditional on those producers and their analytic premises.
"""
import argparse
import hashlib
import json
from pathlib import Path
import zlib


def read_json(path):
    raw = path.read_bytes()
    return json.loads(zlib.decompress(raw) if path.suffix == '.zlib' else raw), hashlib.sha256(raw).hexdigest()


def run(args):
    target, target_hash = read_json(args.target)
    pivots, pivot_hash = read_json(args.pivots)
    reduced, reduced_hash = read_json(args.reduced)
    if (target['schema'], pivots['schema'], reduced['schema']) != (
            'ACTUAL_C27_STRUCTURED_TARGET_V1',
            'ACTUAL_C27_STRUCTURED_JOINT_INTERVAL_V1',
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1'):
        raise ValueError('Require the complete target and its two existing certificates.')
    if (pivots['target_packet_sha256'], reduced['target_packet_sha256'],
            reduced['pivot_result_sha256']) != (target_hash, target_hash, pivot_hash):
        raise ValueError('Target/pivot/Schur packet identity mismatch.')
    metadata = target['metadata']
    if any(reduced[key] != metadata[key] for key in ('N', 'M', 'last_mode', 'basis', 'input_sha256')):
        raise ValueError('Schur metadata do not identify the same complete target.')
    for name, digest in metadata['input_sha256'].items():
        if hashlib.sha256((Path(__file__).resolve().parent / name).read_bytes()).hexdigest() != digest:
            raise ValueError('Actual input pin mismatch: ' + name)
    if pivots['verdict'] != 'inconclusive-pivot':
        raise ValueError('Require a positive prefix followed by an uncertain pivot.')
    prefix = pivots['pivot_modes'][:-1]
    if len(set(prefix)) != len(prefix) or not all(0 <= n <= metadata['N'] for n in prefix):
        raise ValueError('Invalid positive coordinate block.')
    intervals = pivots['pivot_intervals']
    if len(intervals) != len(prefix) + 1 or not all(int(lo) > 0 and int(hi) >= int(lo) for lo, hi in intervals[:-1]):
        raise ValueError('Every reused prefix pivot must be strictly positive.')
    remaining = [n for n in range(metadata['N'] + 1) if n not in set(prefix)]
    if reduced['reused_positive_coordinate_count'] != len(prefix) or reduced['remaining_modes'] != remaining:
        raise ValueError('Reduced block does not complement the certified positive prefix.')
    matrix = reduced['reduced_schur_intervals']
    if len(matrix) != len(remaining) or any(len(row) != len(remaining) for row in matrix):
        raise ValueError('Incomplete reduced matrix.')
    if any(int(lo) > int(hi) for row in matrix for lo, hi in row):
        raise ValueError('Reversed reduced interval.')
    negative = [n for i, n in enumerate(remaining) if int(matrix[i][i][1]) < 0]
    witness = None
    if negative:
        n = negative[0]
        witness = {'remaining_mode': n, 'diagonal_interval': matrix[remaining.index(n)][remaining.index(n)]}
    result = {
        'schema': 'ACTUAL_C27_STRUCTURED_DIAGONAL_CHECK_V1',
        'target_packet_sha256': target_hash, 'pivot_result_sha256': pivot_hash,
        'reduced_result_sha256': reduced_hash, 'grid_bits': reduced['grid_bits'],
        'reused_positive_coordinate_count': len(prefix), 'remaining_modes': remaining,
        'negative_diagonal_modes': negative, 'witness': witness,
        'verdict': 'sufficient-target-not-positive-semidefinite' if negative else 'no-negative-diagonal-certified',
        'conditions': reduced['conditions'] + ['Sound saved dense-Schur interval enclosures for this same target.'],
        'scope': 'A negative reduced diagonal plus the positive prefix rejects this explicit complete sufficient lower target. '
                 'It does not determine the full reduced inertia, actual Schur/Weil sign or every trial on this support. '
                 'No target, solve, spectrum, scalar or moment producer is rerun; no Lean, cofinal positivity or RH result.',
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['verdict'], 'negative diagonal count', len(negative))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('target', 'pivots', 'reduced', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output path.')
    run(args)


if __name__ == '__main__':
    main()
