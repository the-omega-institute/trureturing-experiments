"""Reuse saved cell energies to reject the fixed step-weight certificate.

The reverse weighted triangle is standard. Only exact rational scalar
postprocessing and integer square roots are new; no Fourier work is run.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
PINS = {
    'full-gamma-prime-potential-budget.json':
        '76f2867d013771d78912af8182a8d722c446691f1c3fd1b78c9e229439b4999c',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
}


def rat(q):
    return [str(q.numerator), str(q.denominator)]


def consume(directory, out):
    started = time.monotonic()
    data = []
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted scalar input: '+name)
        data.append(json.loads(raw))
    saved, budget = data
    if (saved['schema'] != 'ACTUAL_C27_PRIME_POTENTIAL_BUDGET_V1'
            or saved['source_sha256'] !=
                '4b7c1e213e3b5ab1b7b54c48e18ec2949b01df691e400dfe35991fd42737b290'
            or saved['packet_sha256']['full-gamma-surplus-budget.json'] !=
                PINS['full-gamma-surplus-budget.json']
            or saved['source_integer_sha256'] != budget['source_integer_sha256']
            or saved['fixed_parameter'] != ['1', '48']
            or saved['grid_bits'] != 640 or len(saved['cell_results']) != 26):
        raise ValueError('Different fixed source, step weight or energy packet.')
    a_lower, previous = F(0), F(0)
    floors = []
    for row in saved['cell_results']:
        left, right, lower = (F(*map(int, row[key])) for key in
                              ('left', 'right', 'potential_lower'))
        lo, hi = (F(int(x), 1 << 640) for x in row['finite_energy'])
        if left != previous or not left < right <= 1 or lower <= 0 or lo > hi or hi < 0:
            raise ValueError('Invalid disjoint cell energy/weight enclosure.')
        a_lower += max(F(0), lo)/lower
        floors.append(lower)
        previous = right
    if previous != 1:
        raise ValueError('Incomplete step-weight partition.')
    rq = budget['exact_rationals']
    far_norm_upper = F(*map(int, rq['far_unweighted_norm_squared_upper']))
    h_upper = F(*map(int, rq['h_upper']))
    floor = min(floors)
    far_weighted_upper = far_norm_upper/floor
    if far_norm_upper < 0 or h_upper <= 0 or F(*map(int,
            saved['exact_rationals']['global_potential_lower'])) != floor:
        raise ValueError('Different complete far source or step majorant.')
    # floor(sqrt(a_lower)) and ceil(sqrt(far_upper)) are exact integers.
    finite_root_lower = isqrt(a_lower.numerator//a_lower.denominator)
    far_root_upper = isqrt(far_weighted_upper.numerator//far_weighted_upper.denominator)
    if far_root_upper**2 < far_weighted_upper:
        far_root_upper += 1
    reverse_lower = F(max(0, finite_root_lower-far_root_upper)**2)
    target = 47*h_upper
    verdict = ('fixed-step-weight-functional-exceeds-shift-budget' if reverse_lower > target
               else 'inconclusive-fixed-step-weight-reverse-triangle')
    result = dict(schema='ACTUAL_C27_PRIME_POTENTIAL_STEP_OBSTRUCTION_V1',
                  packet_sha256=PINS, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  fixed_parameter=['1', '48'], source_integer_sha256=saved['source_integer_sha256'],
                  exact_rationals={key: rat(value) for key, value in dict(
                      finite_step_weight_energy_lower=a_lower,
                      complete_far_step_weight_upper=far_weighted_upper,
                      complete_step_weight_functional_lower=reverse_lower,
                      target_using_h_upper=target).items()},
                  finite_weighted_norm_integer_lower=finite_root_lower,
                  far_weighted_norm_integer_upper=far_root_upper,
                  verdict=verdict, conditions=saved['conditions'],
                  scope='Conditional on sound same-source cells and complete far cap, integral W*|u|^2 exceeds47h '
                      'for this exact26-cell step majorant W. This excludes that sufficient multiplication-weight '
                      'certificate at1/48 even with any far phase satisfying the cap. It does not lower-bound '
                      'integral |u|^2/v, the optimal Pi budget or actual positivity; no Robin/RH/Lean claim.',
                  old_producers_replayed=0, new_Fourier_work=0,
                  new_actions_inverses_solves_Gram_matrices=0,
                  elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('An input changed during consumption.')
    with out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(verdict, 'lower', reverse_lower, 'seconds', result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output; never overwrite accepted data.')
    consume(args.input_dir.resolve(), args.out)
