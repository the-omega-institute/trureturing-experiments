"""Use the missing full-band Bessel projection for the same surplus source.

The source, reflected153-cell weight, complete far cap and two moments
are reused. One standard Fourier multiplication convolution supplies
the projection of W*uJ onto ALL cosine modes0..32768. No old producer
or matrix determinant, inverse or solve is run.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_poly, ctx, fmpq

from diagonal_floor_witness import interval_ball
from full_gamma_prime_potential_budget import ball, rat
from full_gamma_shifted_projection_interval import affine_enclosure, interval
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, FIRST, LAST = 640, 2049, 32768
PINS = {
    'full-gamma-dual-32768-defect.json.zlib':
        '6f320576338562a4e7338906ca66b25ab49e4a1d4460e87e26165cd9579a1154',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
    'full-gamma-prime-potential-budget.json':
        '76f2867d013771d78912af8182a8d722c446691f1c3fd1b78c9e229439b4999c',
    'full-gamma-prime-potential-refined.json':
        '9bad7a02ad1a69ff53f0b6a253cce1b04603d1877ba5da05811cbbd141f9389c',
    'full-gamma-projection-upper.json.zlib':
        'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'full-gamma-shifted-projection-upper.json.zlib':
        '4e024794d7bd7abdf11b270875f5d2fda8295e2d0ca9e2d4117fe4e95e76ad7a',
}


def weight_cosines(cells, last):
    """Integrate a reflected step function, including every rational jump."""
    result = [arb(0) for _ in range(last+1)]
    previous, mean = F(0), F(0)
    for left, right, weight in cells:
        if not previous == left < right <= 1 or weight <= 0:
            raise ValueError('Incomplete positive reflected step partition.')
        mean += (right-left)*weight
        previous = right
    if previous != 1:
        raise ValueError('Step partition does not cover the full interval.')
    result[0] = ball(mean)
    width = 64
    # Endpoint sines vanish exactly at0 and1. Internal coefficient at t_j
    # is (weight_left-weight_right)*sin(pi*l*t_j)/(pi*l).
    for old, new in zip(cells, cells[1:]):
        if old[1] != new[0]:
            raise ValueError('Step endpoints do not agree.')
        point, jump = old[1], ball(old[2]-new[2])
        offsets = [arb.sin_cos_pi_fmpq(fmpq(point.numerator*r, point.denominator))
                   for r in range(width)]
        for start in range(0, last+1, width):
            sine, cosine = arb.sin_cos_pi_fmpq(
                fmpq(point.numerator*start, point.denominator))
            for r in range(max(1-start, 0), min(width, last+1-start)):
                sr, cr = offsets[r]
                result[start+r] += jump*(sine*cr+cosine*sr)
    pi = arb.pi()
    for k in range(1, last+1):
        result[k] /= k*pi
    return result


def multiply_project(values, first, cosines, last):
    """Toeplitz plus Hankel action for positive input modes first>=1."""
    n = first+len(values)-1
    if first < 1 or not values or not 0 <= last <= n or len(cosines) != 2*n+1:
        raise ValueError('Insufficient Fourier coefficients for the projection.')
    symmetric_source = [arb(0) for _ in range(2*n+1)]
    for k, value in enumerate(values, first):
        symmetric_source[n-k] = value
        symmetric_source[n+k] = value
    symmetric_weight = list(reversed(cosines[1:]))+cosines
    product = arb_poly(symmetric_source)*arb_poly(symmetric_weight)
    center = 3*n
    # Nonconstant normalized cosine coefficients are the Laurent
    # coefficients. The constant basis vector has the extra1/sqrt2.
    return [product[center]/arb(2).sqrt()]+[product[center+k] for k in range(1, last+1)]


def coefficient_export(values):
    result = []
    for value in values:
        if not value.is_finite():
            raise ArithmeticError('Nonfinite coefficient; stop without retry.')
        pair = envelope(value, BITS)
        if pair is None or len(pair) != 2 or int(pair[0]) > int(pair[1]):
            raise ArithmeticError('Invalid coefficient envelope; stop without retry.')
        result.append(pair)
    return result


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
            raise ValueError('Different accepted packet: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    defect, budget, old, refined, projection, shifted = packets.values()
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_FINITE_DUAL_DEFECT_V1', 'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1',
            'ACTUAL_C27_PRIME_POTENTIAL_BUDGET_V1',
            'ACTUAL_C27_REFINED_PRIME_POTENTIAL_BUDGET_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
            'ACTUAL_C27_SHIFTED_PROJECTION_UPPER_V1']
            or any(p['grid_bits'] != BITS for p in (defect, old, refined, projection, shifted))
            or (defect['dual_first'], defect['dual_last']) != (FIRST, LAST)
            or (refined['finite_first'], refined['finite_last']) != (FIRST, LAST)
            or len(old['cell_results']) != 26 or len(refined['refined_cells']) != 128
            or refined['fixed_last_cell_subdivisions'] != 128
            or budget['source_bits'] != 384 or budget['source_row'] != 0
            or budget['source_integer_sha256'] !=
                '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
            or any(p['source_integer_sha256'] != budget['source_integer_sha256']
                   for p in (old, refined, shifted))
            or projection['exact_source_integer_sha256'][0] != budget['source_integer_sha256']
            or projection['remaining_modes'] != shifted['remaining_modes']
            or projection['determinant_grid_bits'] != 11584
            or shifted['determinant_grid_bits'] != 11584
            or shifted['fixed_parameter'] != ['1', '48']):
        raise ValueError('Different source, weight, normalization or ceiling roles.')
    for packet in (budget, old, refined, shifted):
        for name, digest in packet['packet_sha256'].items():
            if name in PINS and PINS[name] != digest:
                raise ValueError('Different inherited packet correspondence: '+name)
    if defect['metadata']['basis'] != '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase':
        raise ValueError('Different cosine phase or normalization.')
    cells = [(F(*map(int, row['left'])), F(*map(int, row['right'])),
              1/F(*map(int, row['potential_lower']))) for row in old['cell_results'][:-1]]
    cells += [(F(*map(int, row['left'])), F(*map(int, row['right'])),
               F(*map(int, row['weight']))) for row in refined['refined_cells']]
    maximum = max(weight for _, _, weight in cells)
    if maximum != F(*map(int, refined['exact_rationals']['maximum_weight'])):
        raise ValueError('Different saved complete maximum weight.')
    profile = defect['effective_dual_action']
    if len(profile) != LAST-FIRST+1 or any(len(row) != 1 for row in profile):
        raise ValueError('Incomplete accepted source profile.')
    values = [interval_ball(row[0], BITS) for row in profile]
    a1 = interval_ball(refined['scalar_bounds']['finite_weighted_first_moment'], BITS)
    a2 = interval_ball(refined['scalar_bounds']['finite_weighted_second_moment'], BITS)
    if not a1 > 0 or not a2 > 0:
        raise ValueError('Require the accepted positive same-weight moments.')
    rq = budget['exact_rationals']
    far = F(*map(int, rq['far_unweighted_norm_squared_upper']))
    hlo, hhi = (F(*map(int, rq[key])) for key in ('h_lower', 'h_upper'))
    if far < 0 or not 0 < hlo <= hhi:
        raise ValueError('Invalid complete far cap or h enclosure.')
    print('new step Fourier coefficients0..65536', flush=True)
    cosines = weight_cosines(cells, 2*LAST)
    projected = multiply_project(values, FIRST, cosines, LAST)
    projected_energy = sum((value*value for value in projected), arb(0))
    tail = a2.upper()-max(arb(0), projected_energy.lower())
    if not tail >= 0:
        raise ArithmeticError('Inconsistent Bessel projection and saved a2; stop.')
    cross = 2*ball(far).sqrt()*tail.sqrt()
    complete = a1.upper()+cross.upper()+ball(maximum*far)
    pair = envelope(complete, BITS)
    if pair is None:
        raise ArithmeticError('Nonfinite complete upper; stop without retry.')
    upper = F(int(pair[1]), 1 << BITS)
    if upper <= 0:
        raise ArithmeticError('Nonpositive complete upper.')
    cap = hlo/(hlo+upper)
    d0, d1 = (interval(p['determinant_interval'], 11584) for p in (projection, shifted))
    slope = (48*(d1[0]-d0[1]), 48*(d1[1]-d0[0]))
    if not d0[1] < 0 or not 0 < slope[0] <= slope[1]:
        raise ArithmeticError('Archived affine ceiling has no certified positive slope.')
    root = (-d0[1]/slope[1], -d0[0]/slope[0])
    if not 0 < root[0] <= root[1] < 1:
        raise ArithmeticError('Ceiling repair threshold not certified inside(0,1).')
    target = (hlo*(1-root[1])/root[1], hhi*(1-root[0])/root[0])
    determinant = affine_enclosure(cap, d0, d1)
    scalars = dict(projected_weighted_source_norm_squared=projected_energy,
                   high_band_weighted_source_norm_squared_upper=tail,
                   complete_cross_upper=cross, complete_weighted_budget_upper=complete)
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(value is None for value in exported.values()):
        raise ArithmeticError('Nonfinite scalar export; no precision retry.')
    meets_target = upper < target[0]
    verdict = ('conditional-band-budget-below-projection-capacity-target' if meets_target
               else 'inconclusive-band-budget-projection-capacity-target')
    determinant_sign = 1 if determinant[0] > 0 else (-1 if determinant[1] < 0 else 0)
    if meets_target and determinant_sign != 1:
        raise ArithmeticError('Budget target passed but affine diagnostic not positive; stop.')
    cosine_export = coefficient_export(cosines)
    projection_export = coefficient_export(projected)
    result = dict(schema='ACTUAL_C27_BAND_SURPLUS_BUDGET_V1', packet_sha256=PINS,
                  source_sha256=guards[source],
                  invocation_source_sha256={p.name: digest for p, digest in guards.items()},
                  source_integer_sha256=budget['source_integer_sha256'],
                  source_bits=384, source_row=0, working_bits=1024, grid_bits=BITS,
                  finite_first=FIRST, finite_last=LAST, projection_first=0, projection_last=LAST,
                  reflected_weight_cells=len(cells), scalar_bounds=exported,
                  weight_cosine_coefficients=cosine_export,
                  weighted_source_projection=projection_export,
                  complete_budget_below_conservative_capacity_target=meets_target,
                  ceiling_determinant_at_cap_sign=determinant_sign,
                  exact_rationals=dict(complete_far_norm_squared_upper=rat(far),
                      maximum_weight=rat(maximum), supplied_parameter_cap=rat(cap),
                      ceiling_zero_lower=rat(root[0]), ceiling_zero_upper=rat(root[1]),
                      ceiling_capacity_target_lower=rat(target[0]),
                      ceiling_capacity_target_upper=rat(target[1]),
                      determinant_at_cap_lower=rat(determinant[0]),
                      determinant_at_cap_upper=rat(determinant[1])),
                  verdict=verdict, conditions=refined['conditions']+[
                      'The SAME complete uF is orthogonal to EVERY unphased cosine mode0..32768, not only uJ.',
                      'The saved a2 is the full squared norm of the SAME step-weight product W*uJ; Bessel subtracts a LOWER projection energy.',
                      'All inherited actual common forms, source, legal E, short, unitary and entire-far assumptions apply.',
                      'Affine ceiling propagation is for the SAME exact U2+tau*rV*rV/h; positive upper determinant does not certify a lower or actual positivity.'],
                  scope='One fixed full-band multiplication projection, using standard Fourier '
                      'convolution/Bessel/Cauchy. A better conditional budget supplies eligibility '
                      'only. It gives no sufficient actual score, retained/cofinal positivity, '
                      'Robin, RH, Lean or generic novelty result.',
                  new_weight_cosine_coefficients=2*LAST+1,
                  new_weighted_source_projection_coefficients=LAST+1,
                  new_polynomial_products=1, old_producers_replayed=0,
                  new_determinant_calls=0, parameter_source_precision_scan=False,
                  elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input changed: '+name)
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Invoked source/helper changed: '+path.name)
    with out.open('xb') as stream:
        stream.write(zlib.compress(json.dumps(result, separators=(',', ':')).encode(), 9))
    for key, value in scalars.items():
        print(key, value.str(18), flush=True)
    print('supplied parameter cap', ball(cap).str(18), flush=True)
    print('ceiling determinant at cap', ball(determinant[0]).str(18), ball(determinant[1]).str(18), flush=True)
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use fresh .zlib output; accepted archives are never overwritten.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
