"""Extend the two saved determinant enclosures across the supplied budget.

A rank-one determinant is affine in the shift. Only exact rational
interpolation and the already supplied resolvent budget cap are new;
no matrix determinant, solve, source, action or other producer is run.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from full_gamma_shifted_projection_upper import PINS as ORIGINAL_PINS

HERE = Path(__file__).resolve().parent
PINS = dict(ORIGINAL_PINS, **{
    'full-gamma-shifted-projection-upper.json.zlib':
        '4e024794d7bd7abdf11b270875f5d2fda8295e2d0ca9e2d4117fe4e95e76ad7a',
})
TAU0 = F(1, 48)


def rational(value):
    return [str(value.numerator), str(value.denominator)]


def interval(pair, bits):
    if len(pair) != 2 or int(pair[0]) > int(pair[1]):
        raise ValueError('Invalid archived endpoints.')
    return tuple(F(int(x), 1 << bits) for x in pair)


def affine_enclosure(t, original, shifted):
    """Keep the reversed endpoint direction when extrapolating past tau0."""
    ratio = t/TAU0
    terms = [(1-ratio, original), (ratio, shifted)]
    return tuple(sum((min(c*p[0], c*p[1]) if end == 0 else
                      max(c*p[0], c*p[1]) for c, p in terms), F(0))
                 for end in (0, 1))


def display(value):
    with localcontext() as context:
        context.prec = 24
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    guards = {source: hashlib.sha256(source.read_bytes()).hexdigest()}
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            path = Path(filename).resolve()
            guards[path] = hashlib.sha256(path.read_bytes()).hexdigest()
    packets = []
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted packet: '+name)
        packets.append(json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw))
    old, budget, eligibility, shifted = packets
    if (old['schema'] != 'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1'
            or budget['schema'] != 'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1'
            or eligibility['schema'] != 'ACTUAL_C27_REFINED_RESOLVENT_BUDGET_V1'
            or shifted['schema'] != 'ACTUAL_C27_SHIFTED_PROJECTION_UPPER_V1'
            or old['grid_bits'] != 640 or shifted['grid_bits'] != 640
            or eligibility['grid_bits'] != 640
            or any(p['fixed_parameter'] != ['1', '48'] for p in (eligibility, shifted))
            or eligibility['verdict'] != 'conditional-fixed-shift-budget-admitted'
            or shifted['verdict'] != 'conditional-shifted-comparison-family-excluded'
            or old['determinant_grid_bits'] != 11584
            or shifted['determinant_grid_bits'] != 11584
            or old['remaining_modes'] != shifted['remaining_modes']):
        raise ValueError('Different determinant/shift/budget roles or scales.')
    for name, digest in ORIGINAL_PINS.items():
        if shifted['packet_sha256'][name] != digest:
            raise ValueError('Different saved fixed-shift inputs.')
    if (budget['packet_sha256']['full-gamma-projection-upper.json.zlib'] !=
            PINS['full-gamma-projection-upper.json.zlib']
            or eligibility['packet_sha256']['full-gamma-surplus-budget.json'] !=
                PINS['full-gamma-surplus-budget.json']
            or old['exact_source_integer_sha256'][0] != budget['source_integer_sha256']
            or any(p['source_integer_sha256'] != budget['source_integer_sha256']
                   for p in (eligibility, shifted))
            or shifted['source_bits'] != 384 or shifted['source_row'] != 0):
        raise ValueError('Different common source or inverse-form correspondence.')
    hlo = F(*map(int, budget['exact_rationals']['h_lower']))
    upper = interval(eligibility['scalar_bounds']['complete_weighted_budget_upper'], 640)[1]
    if not hlo > 0 or not upper > 0:
        raise ValueError('Require positive h lower and complete budget upper.')
    cap = hlo/(hlo+upper)
    if not TAU0 < cap < 1:
        raise ArithmeticError('The supplied refined budget does not extend past tau0; stop.')
    d0 = interval(old['determinant_interval'], 11584)
    d1 = interval(shifted['determinant_interval'], 11584)
    endpoint = affine_enclosure(cap, d0, d1)
    ceiling = max(d0[1], endpoint[1])
    result = dict(schema='ACTUAL_C27_SHIFTED_PROJECTION_INTERVAL_V1', packet_sha256=PINS,
                  source_sha256=guards[source],
                  invocation_source_sha256={p.name: h for p, h in guards.items()},
                  source_integer_sha256=budget['source_integer_sha256'],
                  source_bits=384, source_row=0,
                  fixed_reference_parameter=rational(TAU0),
                  exact_rationals=dict(h_lower=rational(hlo), complete_budget_upper=rational(upper),
                      supplied_parameter_cap=rational(cap),
                      determinant_at_cap_lower=rational(endpoint[0]),
                      determinant_at_cap_upper=rational(endpoint[1]),
                      uniform_determinant_upper=rational(ceiling)),
                  verdict=('conditional-supplied-shift-interval-comparison-excluded' if ceiling < 0
                           else 'inconclusive-supplied-shift-interval-ceiling'),
                  conditions=shifted['conditions']+[
                      'The saved complete resolvent upper bounds the SAME p_v independently of tau; it and h_lower supply eligibility for every 0<tau<=h_lower/(h_lower+upper).',
                      'The SAME exact matrix function U2+tau*rV*rV/h has determinant affine in tau; the archived zero and1/48 determinant intervals enclose its two values.',
                      'Affine interpolation/extrapolation uses directed endpoints for both signs of each coefficient; a strict negative endpoint maximum controls the entire interval.'],
                  scope='Only the supplied sufficient eligibility interval for this fixed '
                      'source and M_tau comparison. This is not the maximal actual capacity '
                      'or an exclusion for every actually eligible shift. No actual Weil '
                      'negativity, retained/cofinal positivity, Robin, RH, Lean or novelty claim.',
                  new_determinant_inverse_solve_calls=0, old_producers_replayed=0,
                  parameter_source_precision_scan=False, elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input changed: '+name)
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Invoked source/helper changed: '+path.name)
    with out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('supplied parameter cap', display(cap), flush=True)
    print('uniform determinant upper', display(ceiling), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.json':
        parser.error('Use fresh .json output; accepted archives are never overwritten.')
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
