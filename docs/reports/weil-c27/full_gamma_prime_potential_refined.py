"""Refine one saved potential cell for the same full source at shift 1/48.

The first 25 cell masses and the complete far cap are reused. Only 127
missing interior primitives and the new common-weight moments are formed.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
import time
import zlib

from flint import arb, ctx

from diagonal_floor_witness import interval_ball
from full_gamma_prime_potential_budget import (
    B, BITS, FIRST, LAST, ROW, ball, cells_from_cache, cosine_square,
    rat, sine_primitive,
)
from full_gamma_surplus_shift_gate import minimum
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
SUBCELLS = 128
PINS = {
    'high-floor.json': '3e85480a54d57e2e842964d9be764c755aedfeefd3ff6feec3f851dc0acf8120',
    'full-gamma-dual-32768-defect.json.zlib':
        '6f320576338562a4e7338906ca66b25ab49e4a1d4460e87e26165cd9579a1154',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
    'full-gamma-prime-potential-budget.json':
        '76f2867d013771d78912af8182a8d722c446691f1c3fd1b78c9e229439b4999c',
}


def quadratic_on(cache, row):
    """Restore the accepted rational lower polynomial; no log/Gamma supplier."""
    grid = cache['rational_grid_denominator']
    logs = {int(p): tuple(map(F, pair)) for p, pair in cache['logarithm_intervals'].items()}
    weights, shifts = {}, {}
    for n in row['negative']+row['positive']:
        p, j = cache['prime_powers'][str(n)]
        if n != p**j or not 2 <= n < 27:
            raise ValueError('Different accepted prime-power labels.')
        root = isqrt(n*grid*grid)
        weights[n] = logs[p][1]/F(root, grid)
        shifts[n] = ((F(2*j, 3),)*2 if p == 3 else
                     (2*j*logs[p][0]/(3*logs[3][1]),
                      2*j*logs[p][1]/(3*logs[3][0])))
    ww = sum((weights[n] for n in row['negative']+row['positive']), F(0))
    bb = (sum((weights[n]*shifts[n][0] for n in row['negative']), F(0))
          - sum((weights[n]*shifts[n][1] for n in row['positive']), F(0)))
    dd = sum((weights[n]*shifts[n][1]**2 for n in row['negative']+row['positive']), F(0))
    return ROW-ww-B*dd, 2*B*bb, B*(ROW-ww)


def minimum_on(polynomial, left, right):
    c0, c1, c2 = polynomial
    length = right-left
    _, value = minimum([c0+c1*left+c2*left*left,
                        length*(c1+2*c2*left), c2*length*length])
    return value


def nonnegative_energy(value):
    if not value.is_finite() or value.upper() < 0:
        raise ArithmeticError('Nonfinite or inconsistent original cell mass.')
    pair = envelope(value, BITS)
    if pair is None or int(pair[1]) < 0:
        raise ArithmeticError('Nonfinite or inconsistent nonnegative cell mass.')
    # Intersect only with the independently known mathematical nonnegativity.
    pair[0] = str(max(0, int(pair[0])))
    return pair, interval_ball(pair, BITS)


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    # Invocation provenance only; code bytes are not data-cache reuse keys.
    source_guards = {source: hashlib.sha256(source.read_bytes()).hexdigest()}
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            path = Path(filename).resolve()
            source_guards[path] = hashlib.sha256(path.read_bytes()).hexdigest()
    packets = {}
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted data: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    cache, defect, budget, energy = packets.values()
    if ((cache['cutoff'], cache['weight'], cache['row_upper']) !=
            (27, '1+(29/50)(u/a)^2', '59/10')
            or defect['schema'] != 'ACTUAL_C27_FINITE_DUAL_DEFECT_V1'
            or budget['schema'] != 'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1'
            or energy['schema'] != 'ACTUAL_C27_PRIME_POTENTIAL_BUDGET_V1'
            or energy['packet_sha256'] != {key: PINS[key] for key in list(PINS)[:3]}
            or defect['grid_bits'] != BITS or energy['grid_bits'] != BITS
            or (defect['dual_first'], defect['dual_last']) != (FIRST, LAST)
            or (energy['finite_first'], energy['finite_last']) != (FIRST, LAST)
            or energy['fixed_parameter'] != ['1', '48']
            or budget['source_bits'] != 384 or budget['source_row'] != 0
            or budget['source_integer_sha256'] !=
                '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
            or energy['source_integer_sha256'] != budget['source_integer_sha256']
            or budget['packet_sha256']['full-gamma-dual-32768-defect.json.zlib'] !=
                PINS['full-gamma-dual-32768-defect.json.zlib']
            or defect['metadata']['basis'] !=
                '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'):
        raise ValueError('Different same-source, phase, weight or cell correspondence.')
    old = cells_from_cache(cache)
    saved = energy['cell_results']
    if len(saved) != len(old):
        raise ValueError('Missing accepted cell masses.')
    for row, (left, right, lower) in zip(saved, old):
        if (F(*map(int, row['left'])), F(*map(int, row['right'])),
                F(*map(int, row['potential_lower']))) != (left, right, lower):
            raise ValueError('Different accepted cell correspondence.')
    a0 = interval_ball(energy['scalar_bounds']['finite_unweighted_norm_squared'], BITS)
    if not a0 > 0:
        raise ValueError('Require the saved strictly positive finite norm.')
    rq = budget['exact_rationals']
    far = F(*map(int, rq['far_unweighted_norm_squared_upper']))
    hlow, hupp = (F(*map(int, rq[key])) for key in ('h_lower', 'h_upper'))
    if not 0 < hlow <= hupp or far < 0:
        raise ValueError('Require positive h and the complete nonnegative far cap.')
    target = 47*hlow
    left = old[-1][0]
    polynomial = quadratic_on(cache, cache['support_cells'][-1])
    if minimum_on(polynomial, left, F(1)) != F(cache['support_cells'][-1]['polynomial_minimum']):
        raise ValueError('Restored polynomial does not match the archived minimum.')
    points = [left+(1-left)*F(j, SUBCELLS) for j in range(SUBCELLS+1)]
    refined = []
    for low, high in zip(points, points[1:]):
        margin = minimum_on(polynomial, low, high)
        if margin <= 0:
            raise ArithmeticError('Nonpositive restored subcell margin.')
        weight = (1+B*high*high)/margin
        if weight > 1/old[-1][2]:
            raise ArithmeticError('The refined majorant increased above the old weight.')
        refined.append((low, high, weight))
    weights = [1/lower for _, _, lower in old[:-1]]+[w for _, _, w in refined]
    maximum = max(weights)
    far_weighted = maximum*far
    result_cells, scalars = [], dict(finite_unweighted_norm_squared=a0)
    missing_primitives = 0
    if far_weighted >= target:
        verdict = 'fixed-refined-recipe-stopped-at-complete-far-cap'
    else:
        profile = defect['effective_dual_action']
        if len(profile) != LAST-FIRST+1 or any(len(row) != 1 for row in profile):
            raise ValueError('Incomplete original finite source profile.')
        values = [interval_ball(row[0], BITS) for row in profile]
        # This harmonic intermediate was not archived. Restore it with the
        # existing algorithm; reuse a0 rather than contracting another norm.
        harmonic = cosine_square(values, FIRST)
        pi = arb.pi()
        primitive_coefficients = [arb(0)]+[v/(k*pi) for k, v in enumerate(harmonic[1:], 1)]
        endpoints = {left: a0*ball(1-left)-interval_ball(saved[-1]['finite_energy'], BITS),
                     F(1): arb(0)}
        for point in points[1:-1]:
            endpoints[point] = sine_primitive(primitive_coefficients, point)
            missing_primitives += 1
        a1, a2, refined_mass = arb(0), arb(0), arb(0)
        for row, weight in zip(saved[:-1], weights[:25]):
            _, mass = nonnegative_energy(interval_ball(row['finite_energy'], BITS))
            a1 += ball(weight)*mass
            a2 += ball(weight*weight)*mass
        for low, high, weight in refined:
            pair, mass = nonnegative_energy(a0*ball(high-low)+endpoints[high]-endpoints[low])
            refined_mass += mass
            a1 += ball(weight)*mass
            a2 += ball(weight*weight)*mass
            result_cells.append(dict(left=rat(low), right=rat(high), weight=rat(weight),
                                     finite_energy=pair))
        if not refined_mass.overlaps(interval_ball(saved[-1]['finite_energy'], BITS)):
            raise ArithmeticError('Refined masses are inconsistent with the accepted total.')
        if not a1.lower() >= 0:
            raise ArithmeticError('Nonnegative first moment enclosure required.')
        variance_upper = a2.upper()-a1.lower()**2/a0.upper()
        if not variance_upper >= 0:
            raise ArithmeticError('Inconsistent centered-variance upper.')
        cross_upper = 2*ball(far).sqrt()*variance_upper.sqrt()
        complete = a1.upper()+cross_upper+ball(far_weighted)
        scalars.update(finite_weighted_first_moment=a1, finite_weighted_second_moment=a2,
                       refined_last_cell_mass=refined_mass, centered_variance_upper=variance_upper,
                       complete_cross_upper=cross_upper, complete_weighted_budget_upper=complete)
        verdict = ('conditional-fixed-shift-budget-admitted' if complete.upper() <= ball(target)
                   else 'inconclusive-fixed-refined-centered-budget')
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(value is None for value in exported.values()):
        raise ArithmeticError('Nonfinite scalar; no precision/cell retry.')
    result = dict(schema='ACTUAL_C27_REFINED_PRIME_POTENTIAL_BUDGET_V1', packet_sha256=PINS,
                  source_sha256=source_guards[source],
                  invocation_source_sha256={path.name: digest for path, digest in source_guards.items()},
                  working_bits=1024, grid_bits=BITS, fixed_parameter=['1', '48'],
                  source_integer_sha256=budget['source_integer_sha256'],
                  finite_first=FIRST, finite_last=LAST, reused_unchanged_cells=25,
                  fixed_last_cell_subdivisions=SUBCELLS, new_interior_primitives=missing_primitives,
                  restored_unarchived_harmonic_intermediates=int(missing_primitives > 0),
                  restored_last_cell_polynomial=[rat(x) for x in polynomial],
                  exact_rationals=dict(maximum_weight=rat(maximum),
                       complete_far_weighted_upper=rat(far_weighted),
                       target_using_h_lower=rat(target)),
                  scalar_bounds=exported, refined_cells=result_cells, verdict=verdict,
                  conditions=energy['conditions']+[
                      'The restored last-cell polynomial is the same accepted lower for q*v on the complete outward cell.',
                      'The original full source splits into unweighted orthogonal uJ and uF; all new moments use one common reflected majorant.',
                  ],
                  scope='Conditional same-source fixed1/48 budget application only. Failure of this '
                       'chosen complete upper is inconclusive, not actual ineligibility or a Weil sign. '
                       'Actual sufficient score, retained/cofinal positivity, Robin, RH and Lean remain open.',
                  old_producers_replayed=0, old_cell_integrals_recomputed=0,
                  new_actions_inverses_solves_Grams=0, parameter_scan=False,
                  cell_or_precision_retry=False, elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input changed during consumption.')
    for path, digest in source_guards.items():
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
