"""Apply cached weighted rows to the same full source at the fixed shift1/48.

FLINT polynomial multiplication and anchored trigonometric sums supply
only the missing finite weighted source energy. Every unknown far cross
term is paid by the weighted triangle inequality. No old producer is run.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time
import zlib

from flint import arb, arb_poly, ctx, fmpq

from diagonal_floor_witness import interval_ball
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
PINS = {
    'high-floor.json': '3e85480a54d57e2e842964d9be764c755aedfeefd3ff6feec3f851dc0acf8120',
    'full-gamma-dual-32768-defect.json.zlib':
        '6f320576338562a4e7338906ca66b25ab49e4a1d4460e87e26165cd9579a1154',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
}
FIRST, LAST, BITS = 2049, 32768, 640
B, ROW = F(29, 50), F(59, 10)


def ball(q):
    q = F(q)
    return arb(fmpq(q.numerator, q.denominator))


def rat(q):
    return [str(q.numerator), str(q.denominator)]


def cells_from_cache(cache):
    cells = []
    previous = F(0)
    for row in cache['support_cells']:
        left, right = F(row['left']), F(row['right'])
        margin = F(row['polynomial_minimum'])
        if left != previous or not 0 <= left < right <= 1 or margin <= F(1, 400):
            raise ValueError('Different accepted outward cell partition or margin.')
        # The accepted lower is q*v >= margin on this entire rational cell.
        lower = margin/(1+B*right*right)
        cells.append((left, right, lower))
        previous = right
    if previous != 1 or len(cells) != 26:
        raise ValueError('Incomplete accepted26-cell cover.')
    return cells


def cosine_square(coefficients, first):
    """2*(sum a_n cos(pi*n*x))^2, in its cosine basis."""
    count = len(coefficients)
    poly = arb_poly(coefficients)
    correlation = poly*arb_poly(list(reversed(coefficients)))
    square = poly*poly
    result = [correlation[count-1]]+[arb(0)]*(2*(first+count-1))
    for k in range(1, count):
        result[k] += 2*correlation[count-1+k]
    for k in range(2*count-1):
        result[2*first+k] += square[k]
    return result


def sine_primitive(harmonic_over_k_pi, point):
    """Use exact rational anchors and finite trigonometric addition blocks."""
    width = 64
    offsets = [arb.sin_cos_pi_fmpq(fmpq(point.numerator*r, point.denominator))
               for r in range(width)]
    total = arb(0)
    for first in range(0, len(harmonic_over_k_pi), width):
        sine, cosine = arb.sin_cos_pi_fmpq(fmpq(point.numerator*first, point.denominator))
        for r in range(min(width, len(harmonic_over_k_pi)-first)):
            sr, cr = offsets[r]
            total += harmonic_over_k_pi[first+r]*(sine*cr+cosine*sr)
    return total


def weighted_finite_energy(coefficients, cells, norm):
    harmonic = cosine_square(coefficients, FIRST)
    pi = arb.pi()
    primitive_coefficients = [arb(0)]+[value/(k*pi) for k, value in enumerate(harmonic[1:], 1)]
    endpoints = {F(0): arb(0), F(1): arb(0)}
    for point in sorted({x for left, right, _ in cells for x in (left, right)}):
        if point not in endpoints:
            endpoints[point] = sine_primitive(primitive_coefficients, point)
    total, exported = arb(0), []
    for left, right, lower in cells:
        energy = norm*ball(right-left)+endpoints[right]-endpoints[left]
        if not energy.upper() >= 0:
            raise ArithmeticError('A nonnegative cell energy has no nonnegative enclosure.')
        # Interval dependence can cross zero; only its certified upper is used.
        total += energy.upper()/ball(lower)
        exported.append(dict(left=rat(left), right=rat(right), potential_lower=rat(lower),
                             finite_energy=envelope(energy, BITS)))
    return total.upper(), exported


def consume(directory, out):
    started = time.monotonic()
    packets = {}
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    cache, defect, budget = packets.values()
    if (cache['cutoff'], cache['weight'], cache['row_upper']) != (27, '1+(29/50)(u/a)^2', '59/10'):
        raise ValueError('Different accepted weighted row.')
    if (defect['schema'] != 'ACTUAL_C27_FINITE_DUAL_DEFECT_V1'
            or budget['schema'] != 'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1'
            or defect['grid_bits'] != BITS
            or (defect['dual_first'], defect['dual_last']) != (FIRST, LAST)
            or budget['packet_sha256']['full-gamma-dual-32768-defect.json.zlib'] !=
                PINS['full-gamma-dual-32768-defect.json.zlib']
            or budget['source_integer_sha256'] !=
                '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
            or budget['source_bits'] != 384 or budget['source_row'] != 0
            or defect['metadata']['basis'] !=
                '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'):
        raise ValueError('Different original full unscaled action or cosine convention.')
    profile = defect['effective_dual_action']
    if len(profile) != LAST-FIRST+1 or any(len(row) != 1 for row in profile):
        raise ValueError('Incomplete finite source profile.')
    cells = cells_from_cache(cache)
    floor = min(lower for _, _, lower in cells)
    if floor < F(1, 632):
        raise ValueError('The cached cell floor must retain the inherited global floor.')
    rq = budget['exact_rationals']
    far = F(*map(int, rq['far_unweighted_norm_squared_upper']))
    hlow, hupp = (F(*map(int, rq[key])) for key in ('h_lower', 'h_upper'))
    if not 0 < hlow <= hupp or far < 0:
        raise ValueError('Require accepted positive h and complete far cap.')
    values = [interval_ball(row[0], BITS) for row in profile]
    norm = sum((value**2 for value in values), arb(0))
    if not norm > 0:
        raise ArithmeticError('Require nonzero finite source norm.')
    target = 47*hlow
    pfar = far/floor
    minimum_recipe = (ball(pfar).sqrt()+(norm/ball(ROW)).lower().sqrt())**2
    scalars = dict(finite_unweighted_norm_squared=norm,
                   cheap_finite_weighted_energy_lower=(norm/ball(ROW)).lower(),
                   cheap_triangle_recipe_lower=minimum_recipe)
    if pfar < target:
        scalars['finite_weighted_energy_allowance'] = (ball(target).sqrt()-ball(pfar).sqrt())**2
    cell_results = []
    if pfar >= target or minimum_recipe > ball(target):
        verdict = 'fixed-triangle-recipe-stopped-before-weighted-energy'
    else:
        pj, cell_results = weighted_finite_energy(values, cells, norm)
        complete = (pj.sqrt()+ball(pfar).sqrt())**2
        scalars.update(finite_weighted_energy_upper=pj, complete_weighted_budget_upper=complete)
        verdict = ('conditional-fixed-shift-budget-admitted' if complete <= ball(target)
                   else 'inconclusive-fixed-triangle-weighted-budget')
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(value is None for value in exported.values()):
        raise ArithmeticError('Nonfinite scalar; stop without a precision/cell retry.')
    result = dict(schema='ACTUAL_C27_PRIME_POTENTIAL_BUDGET_V1', packet_sha256=PINS,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  working_bits=1024, grid_bits=BITS, fixed_parameter=['1', '48'],
                  finite_first=FIRST, finite_last=LAST, source_integer_sha256=budget['source_integer_sha256'],
                  exact_rationals=dict(global_potential_lower=rat(floor),
                       complete_far_weighted_upper=rat(pfar), target_using_h_lower=rat(target)),
                  scalar_bounds=exported, cell_results=cell_results, verdict=verdict,
                  conditions=budget['conditions']+[
                      'The actual global difference dominates multiplication by v=59/10-r_q on the same common closed form domain.',
                      'The original unphased cosine identification is unitary and E/O orthogonal; finite E shifts preserve the form domain.',
                      'Pi is the same finite-E short; the full u includes the complete unknown far action.',
                      'The sound accepted26-cell polynomial margins lower-bound q*v everywhere up to null endpoints.',
                  ],
                  scope='The fixed source budget only, conditional on the listed actual identifications and sound caches. '
                      'A failed chosen upper or stopped triangle recipe is not actual ineligibility or negativity. '
                      'No actual retained positivity, good-block Gamma upper, cofinal result, Robin tail, RH or Lean claim.',
                  old_producers_replayed=0, new_actions_inverses_solves_Grams=0,
                  parameter_scan=False, cell_or_precision_retry=False,
                  elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('An accepted input changed during consumption.')
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
        parser.error('Use a fresh output; never overwrite accepted data.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
