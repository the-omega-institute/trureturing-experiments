"""Generate c9 scalar inputs using Liu's published, source-only oracle.

The attributed checker is an external input, never vendored. Its complete
6400-bit calculation includes the finite Binet remainder. No author moment
packet, retained matrix, trial inverse or positivity certificate is read.
"""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path


CHECKER_SHA = '4ed7661f1f869ed6d37d7bc52f1a0889a8d305a790663dc0ae43d76f2ea0e356'
CHECKER_URL = ('https://raw.githubusercontent.com/luciferyu666/'
               'certified-weil-positivity/'
               'b6cd2183c1e79c6c27a34267812a7b2d73ed1b59/'
               'frozen-source/reproduction/source/01-check.py')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checker', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('Use a fresh output path.')
    source = args.checker.read_bytes()
    if hashlib.sha256(source).hexdigest() != CHECKER_SHA:
        raise ValueError('Unpinned external moment oracle.')
    namespace = {'__file__': str(args.checker.resolve()),
                 '__name__': 'liu_source_moment_oracle'}
    exec(compile(source, str(args.checker), 'exec'), namespace)
    if namespace['BITS'] != 6400:
        raise ValueError('Unexpected oracle precision.')
    started = time.monotonic()
    print('Computing the complete published 6400-bit moment oracle.', flush=True)
    oracle = namespace['complete_oracle']()
    scale = 1 << (6400 - 1024)
    rows = [{'q': q, 'lo': str(lo // scale), 'hi': str(-(-hi // scale))}
            for q, (lo, hi) in enumerate(oracle)]
    if len(rows) != 1024 or any(not 0 <= int(r['hi']) - int(r['lo']) <= 2
                                for r in rows):
        raise ValueError('Incomplete or excessively wide moment enclosure.')
    namespace['check_rows'](rows, oracle)
    controls = namespace['mutation_controls'](rows, oracle)
    if args.checker.read_bytes() != source:
        raise ValueError('Moment oracle changed during execution.')
    result = {
        'schema': 'C9_VERIFIED_SOURCE_MOMENTS_V1',
        'parameters': {'omega': 256, 'beta': '7/2', 'outputBits': 1024},
        'definition': '(256/pi)*integral_0^1 x^(2q)*(A(256x)-7/2) dx; '
                      'A(t)=Re digamma(1/4+it/2)-log pi',
        'source_checker': {'url': CHECKER_URL, 'sha256': CHECKER_SHA,
                           'working_bits': 6400},
        'source_proof': 'Liu, Appendix D.1, Integrated digamma moments',
        'checked_rows': len(rows),
        'maximum_grid_width': max(int(r['hi']) - int(r['lo']) for r in rows),
        'midpoint_error_upper_power_of_two': -1024,
        'complete_oracle_interval_containment': True,
        'binet_remainder_included_once': True,
        'mutation_controls': controls,
        'python': platform.python_version(),
        'elapsed_seconds': time.monotonic() - started,
        'rows': rows,
        'scope': 'Source-only execution of the attributed existing oracle. '
                 'No independent implementation, matrix sign, Lean or RH proof.'
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}), flush=True)


if __name__ == '__main__':
    main()
