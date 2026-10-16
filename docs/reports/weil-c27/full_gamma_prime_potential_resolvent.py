"""Use the saved refined masses in one fixed orthogonal resolvent bound.

No coefficient, Fourier primitive, source/action or previous moment is
recomputed. The auxiliary Cauchy price is fixed at 640, above the inherited
uniform weight bound 632; the actual shift stays 1/48.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time

from flint import arb, ctx

from diagonal_floor_witness import interval_ball
from full_gamma_prime_potential_budget import BITS, ball, rat
from full_gamma_prime_potential_refined import nonnegative_energy
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
PRICE = F(640)
PINS = {
    'full-gamma-prime-potential-budget.json':
        '76f2867d013771d78912af8182a8d722c446691f1c3fd1b78c9e229439b4999c',
    'full-gamma-prime-potential-refined.json':
        '9bad7a02ad1a69ff53f0b6a253cce1b04603d1877ba5da05811cbbd141f9389c',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
}


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    guards = {source: hashlib.sha256(source.read_bytes()).hexdigest()}
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            path = Path(filename).resolve()
            guards[path] = hashlib.sha256(path.read_bytes()).hexdigest()
    packets = {}
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted data: '+name)
        packets[name] = json.loads(raw)
    old, refined, budget = packets.values()
    if (old['schema'] != 'ACTUAL_C27_PRIME_POTENTIAL_BUDGET_V1'
            or refined['schema'] != 'ACTUAL_C27_REFINED_PRIME_POTENTIAL_BUDGET_V1'
            or budget['schema'] != 'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1'
            or any(p['grid_bits'] != BITS for p in (old, refined))
            or any(p['fixed_parameter'] != ['1', '48'] for p in (old, refined))
            or refined['packet_sha256']['full-gamma-prime-potential-budget.json'] !=
                PINS['full-gamma-prime-potential-budget.json']
            or refined['packet_sha256']['full-gamma-surplus-budget.json'] !=
                PINS['full-gamma-surplus-budget.json']
            or old['packet_sha256']['full-gamma-surplus-budget.json'] !=
                PINS['full-gamma-surplus-budget.json']
            or refined['reused_unchanged_cells'] != 25
            or refined['fixed_last_cell_subdivisions'] != 128
            or refined['new_interior_primitives'] != 127
            or (refined['finite_first'], refined['finite_last']) != (2049, 32768)
            or any(p['source_integer_sha256'] != budget['source_integer_sha256']
                   for p in (old, refined))
            or budget['source_bits'] != 384 or budget['source_row'] != 0
            or budget['source_integer_sha256'] !=
                '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'):
        raise ValueError('Different fixed source or refined weight correspondence.')
    if len(old['cell_results']) != 26 or len(refined['refined_cells']) != 128:
        raise ValueError('Incomplete saved partition.')
    cells = [(row, 1/F(*map(int, row['potential_lower'])))
             for row in old['cell_results'][:-1]]
    cells += [(row, F(*map(int, row['weight']))) for row in refined['refined_cells']]
    maximum = max(weight for _, weight in cells)
    if maximum != F(*map(int, refined['exact_rationals']['maximum_weight'])) or not 0 < maximum < PRICE:
        raise ValueError('Require the same full maximum strictly below fixed price640.')
    previous, resolvent = F(0), arb(0)
    for row, weight in cells:
        left, right = (F(*map(int, row[key])) for key in ('left', 'right'))
        if not previous == left < right <= 1 or weight <= 0:
            raise ValueError('Different or incomplete positive reflected partition.')
        _, mass = nonnegative_energy(interval_ball(row['finite_energy'], BITS))
        resolvent += mass/ball(PRICE-weight)
        previous = right
    if previous != 1 or not resolvent > 0:
        raise ArithmeticError('Incomplete partition or nonpositive resolvent moment.')
    a0 = interval_ball(old['scalar_bounds']['finite_unweighted_norm_squared'], BITS)
    if not a0 > 0:
        raise ArithmeticError('Positive cached finite norm required.')
    rq = budget['exact_rationals']
    far = F(*map(int, rq['far_unweighted_norm_squared_upper']))
    hlow = F(*map(int, rq['h_lower']))
    if far < 0 or hlow <= 0:
        raise ValueError('Require complete nonnegative far cap and positive h.')
    upper = ball(PRICE)*(a0.upper()+ball(far))-a0.lower()**2/resolvent.upper()
    target = 47*hlow
    exported = {key: envelope(value, BITS) for key, value in
                dict(resolvent_moment=resolvent, complete_weighted_budget_upper=upper).items()}
    if any(value is None for value in exported.values()):
        raise ArithmeticError('Nonfinite scalar; stop without any retry.')
    verdict = ('conditional-fixed-shift-budget-admitted' if upper.upper() <= ball(target)
               else 'inconclusive-fixed-price-resolvent-budget')
    result = dict(schema='ACTUAL_C27_REFINED_RESOLVENT_BUDGET_V1', packet_sha256=PINS,
                  source_sha256=guards[source],
                  invocation_source_sha256={p.name: h for p, h in guards.items()},
                  working_bits=1024, grid_bits=BITS, fixed_parameter=['1', '48'],
                  fixed_cauchy_price=rat(PRICE), source_integer_sha256=budget['source_integer_sha256'],
                  reused_cell_masses=len(cells), scalar_bounds=exported,
                  exact_rationals=dict(maximum_weight=rat(maximum),
                      far_unweighted_norm_squared_upper=rat(far), target_using_h_lower=rat(target)),
                  verdict=verdict, conditions=refined['conditions']+[
                      'The SAME full source uJ+uF has unweighted orthogonal finite and complete far parts.',
                      'All saved153 cell masses enclose the same original finite source and reflected refined weight.',
                  ],
                  scope='Standard weighted-Cauchy application at fixed price640 and actual shift1/48; '
                      'a passing upper certifies only the conditional same-source directional budget. '
                      'No actual sufficient score, retained/cofinal positivity, Robin, RH, Lean or novelty claim.',
                  new_coefficient_or_primitive_arithmetic=0, new_actions_inverses_solves_Grams=0,
                  old_producers_replayed=0, price_or_parameter_scan=False,
                  precision_retry=False, elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input changed during consumption.')
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Invoked source/helper changed during consumption: '+path.name)
    with out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output; do not overwrite accepted data.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
