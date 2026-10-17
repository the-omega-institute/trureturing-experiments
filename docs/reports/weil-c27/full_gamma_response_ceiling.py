"""Bound the strength of saved scalar-only lower certificates.

Exact rational endpoint arithmetic on an accepted packet; no scientific
supplier or plane calculation is repeated. The ceiling concerns a
certificate algorithm, never the actual energy or a realized row.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time
import zlib

HERE = Path(__file__).resolve().parent
INPUT = 'full-gamma-cached-source-response.json.zlib'
INPUT_SHA = '2f6fa2f46f91dfa0449d5bec2523cdbd91a0e7587284f8e417ccb5bd2068ee03'


def optimistic_ceiling(a, b, c, lo, hi):
    """Convex p(z)=a*z^2+2*b*z+c is at most its endpoint maximum."""
    if any(type(x) is not Fraction for x in (a, b, c, lo, hi)):
        raise TypeError('Exact Fraction coefficients and endpoints required.')
    if a <= 0 or lo > hi:
        raise ValueError('Positive curvature and an ordered interval required.')
    left, right = a*lo*lo+2*b*lo+c, a*hi*hi+2*b*hi+c
    return left, right, max(left, right)


def read_rational(values, key):
    pair = values[key]
    if (not isinstance(pair, list) or len(pair) != 2
            or any(not isinstance(v, str) for v in pair)):
        raise ValueError('Expected a saved exact-rational pair: '+key)
    numerator, denominator = map(int, pair)
    if denominator <= 0:
        raise ValueError('Positive rational denominator required: '+key)
    return Fraction(numerator, denominator)


def encode(value):
    return [str(value.numerator), str(value.denominator)]


def consume(out):
    started = time.perf_counter()
    source = Path(__file__).resolve()
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    raw = (HERE/INPUT).read_bytes()
    if hashlib.sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Different accepted cached-source response packet.')
    data = json.loads(zlib.decompress(raw))
    if (data['schema'] != 'ACTUAL_C27_CACHED_SOURCE_RESPONSE_V1'
            or (data['source_bits'], data['source_row'], data['grid_bits'],
                data['working_bits'], data['normalization_scale']) != (384, 0, 640, 1024, 32768)
            or data['actual_center_sign_claim'] or data['actual_full_border_sign_claim']
            or data['actual_rh_proof']):
        raise ValueError('Different saved source or certificate scope.')
    values = data['exact_rationals']
    rat = lambda key: read_rational(values, key)
    lo, hi = rat('zeta_lower'), rat('zeta_upper')
    if lo > hi:
        raise ValueError('Reversed saved response interval.')
    recipes = {}
    for recipe in ('fixed', 'energy_endpoint'):
        ends = {}
        for endpoint in ('theta_lo', 'theta_hi'):
            prefix = recipe+'_'+endpoint+'_'
            a = rat(prefix+'quadratic_coefficient_lower')
            b = rat(prefix+'linear_midpoint')
            radius = rat(prefix+'linear_radius_upper')
            c = rat(prefix+'constant_lower')
            if radius < 0:
                raise ValueError('Negative saved linear-radius allowance.')
            left, right, ceiling = optimistic_ceiling(a, b, c, lo, hi)
            ends[endpoint] = dict(
                quadratic_coefficient_lower=encode(a), linear_midpoint=encode(b),
                linear_radius_upper=encode(radius), constant_lower=encode(c),
                optimistic_left_value=encode(left), optimistic_right_value=encode(right),
                optimistic_endpoint_ceiling=encode(ceiling))
        issued_cap = min(read_rational(v, 'optimistic_endpoint_ceiling')
                         for v in ends.values())
        retained = rat(recipe+'_retained_floor_lower')
        retained_cap = max(retained, issued_cap)
        recipes[recipe] = dict(theta_endpoints=ends, current_retained_lower=encode(retained),
            new_scalar_certificate_ceiling=encode(issued_cap),
            retained_scalar_certificate_ceiling=encode(retained_cap),
            scalar_only_refinement_cannot_certify_positive=retained_cap < 0)
    if (hashlib.sha256((HERE/INPUT).read_bytes()).hexdigest() != INPUT_SHA
            or hashlib.sha256(source.read_bytes()).hexdigest() != source_sha):
        raise RuntimeError('Scalar-ceiling inputs changed during this application.')
    result = dict(schema='C27_SAVED_SCALAR_CERTIFICATE_CEILING_V1',
        input_packet=INPUT, input_sha256=INPUT_SHA, source_sha256=source_sha,
        source_integer_sha256=data['source_integer_sha256'],
        source_bits=384, source_row=0, grid_bits=640, working_bits='exact rational',
        response_interval=[encode(lo), encode(hi)], recipes=recipes,
        all_fixed_recipe_scalar_refinements_excluded=all(
            v['scalar_only_refinement_cannot_certify_positive'] for v in recipes.values()),
        new_saved_quadratic_endpoint_evaluations=8, new_recipe_certificate_caps=2,
        new_plane_floor_calls=0, old_producers_or_fixtures_replayed=0,
        new_source_sine_phase_moment_inverse_or_solve_suppliers=0,
        parameter_source_precision_cutoff_scan=False, parameters_changed=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False,
        actual_rh_proof=False, elapsed_seconds=time.perf_counter()-started,
        conditions=[
            'Every nonempty closed J inside the saved response interval is allowed, including a singleton. Original directed coefficients, theta endpoints, k/D, source and form stay fixed.',
            'The original issued quadratic subtracts 2*radius*sup_J|z|. Dropping this nonnegative debit only upper-bounds the issued LOWER certificate, never the actual energy.',
            'Convexity bounds each optimistic quadratic by its endpoint maximum. The minimum over theta endpoints and any J is no larger than the minimum of these maxima.',
            'Retaining earlier certified lowers uses their maximum with the new issued lower. A negative ceiling excludes positive certification only by this unchanged scalar-refinement recipe.',
            'Changed joint constraints, improved coefficient/physical/charge bounds or another projection are outside this diagnostic. No negative actual row or actual paid-energy upper is constructed.',
            'RH, Robin and actual-center/full-border/good-block/cofinal positivity remain unresolved. This is a classical convexity application, not a new generic theorem, novelty or Lean claim.',
        ])
    with Path(out).open('x', encoding='utf-8') as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write('\n')
    print('saved scalar-certificate ceiling supplied; seconds', result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    consume(args.out)
