#!/usr/bin/env python3
"""Periodic floor energy and reciprocal-square samples of existing weights.

Arbitrary clocks A=m**8 are not selected Robin sources. No old response grid
is enumerated. Floating-point energy and Decimal spot checks are numerical
diagnostics, not certified bounds, asymptotics, or a sampling-transfer theorem.
"""

import argparse
from array import array
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

from mixed_sieve_response import (
    decimal_response, float_weights, lcm_coefficients, mobius_table,
)


def primes_up_to(limit):
    composite = bytearray(limit + 1)
    for p in range(2, math.isqrt(limit) + 1):
        if not composite[p]:
            for n in range(p * p, limit + 1, p):
                composite[n] = 1
    return [p for p in range(2, limit + 1) if not composite[p]]


def periodic_energy(h, limit):
    coefficients = array('d', [0.0]) * (limit + 1)
    for d, value in h:
        coefficients[d] = value / d
    jordan = [q * q for q in range(limit + 1)]
    for p in primes_up_to(limit):
        for q in range(p, limit + 1, p):
            jordan[q] -= jordan[q] // (p * p)
    channel_energy = []
    for q in range(1, limit + 1):
        channel = math.fsum(coefficients[d] for d in range(q, limit + 1, q))
        channel_energy.append(jordan[q] * channel * channel / 12)
    variance = math.fsum(channel_energy)
    pair_error = None
    if len(h) <= 100:
        pair_value = math.fsum(
            a * b * math.gcd(d, e) ** 2 / (12 * d * e)
            for d, a in h for e, b in h)
        pair_error = abs(variance - pair_value)
        if pair_error > 1e-11 * max(1.0, variance):
            raise ArithmeticError('Jordan and direct pair energy disagree.')
    center = math.fsum(value for _, value in h) / 2
    return {
        'center': center,
        'centered_variance': variance,
        'centered_standard_deviation': math.sqrt(variance),
        'full_floor_mean_square': variance + center * center,
        'integer_frequency_energy': channel_energy[0],
        'other_rational_frequency_energy': math.fsum(channel_energy[1:]),
        'direct_pair_energy_error': pair_error,
    }


def sample_squares(m, h, mean, mu):
    clock = m ** 8
    cutoff = m ** 2
    samples = []
    smallest = None
    point_count = 0
    for quarter in [4, 5, 6, 7, 8]:
        numerator = quarter * clock
        denominator = 4
        prime_limit = math.isqrt(numerator // (denominator * cutoff))
        costs = []
        for p in primes_up_to(prime_limit):
            if denominator * cutoff * p * p >= numerator:
                continue
            point_count += 1
            y = Fraction(numerator, denominator * p * p)
            floor_correction = math.fsum(
                value * (y.numerator % (y.denominator * d)) /
                (y.denominator * d) for d, value in h)
            response = float(y) * mean - floor_correction
            residual = max(0.0, floor_correction - float(y) * max(mean, 0.0))
            costs.append(math.log(p) * residual)
            if smallest is None or response < smallest['H_float']:
                smallest = {'x': Fraction(numerator, denominator),
                            'p': p, 'y': y, 'H_float': response}
        samples.append({
            'x': str(Fraction(numerator, denominator)),
            'upper_weight_square_residual':
                math.fsum(costs) / (numerator / denominator),
            'prime_count': len(costs),
        })
    if smallest is None:
        raise ArithmeticError('The selected finite scales have no square samples.')
    y = smallest['y']
    _, first = decimal_response(m, mu, [y], 80)
    _, second = decimal_response(m, mu, [y], 120)
    with localcontext() as ctx:
        ctx.prec = 120
        precision_error = abs(first[y] - second[y])
    if precision_error > Decimal('1e-55'):
        raise ArithmeticError('The Decimal spot-check precisions disagree.')
    smallest = {**smallest, 'x': str(smallest['x']), 'y': str(y),
                'H_decimal120': str(second[y]),
                'float_vs_decimal_error': abs(smallest['H_float'] - float(second[y])),
                'decimal80_vs_decimal120_error': str(precision_error)}
    return samples, point_count, smallest


def diagnose(m, existing_means):
    if m not in existing_means:
        raise ValueError('The scale needs an existing high-precision mean record.')
    if m < 2 or m ** 4 > 1_000_000:
        raise ValueError('The scale must be at least 2 with m**4 <= 1,000,000.')
    mu = mobius_table(m ** 3)
    h = lcm_coefficients(*float_weights(m, mu))
    mean = float(existing_means[m])
    rebuilt_mean = math.fsum(value / d for d, value in h)
    if abs(mean - rebuilt_mean) > 1e-12:
        raise ArithmeticError('The reused mean and rebuilt coefficients disagree.')
    energy = periodic_energy(h, m ** 4)
    absolute_sum = math.fsum(abs(value) for _, value in h)
    samples, point_count, smallest = sample_squares(m, h, mean, mu)
    return {
        'm': m, 'A': m ** 8, 'support_upper_bound': m ** 4,
        'grouped_lcm_count': len(h),
        'reused_mean_decimal120': existing_means[m],
        'sum_absolute_grouped_coefficients': absolute_sum,
        'periodic_energy': energy,
        'standard_deviation_over_absolute_sum':
            energy['centered_standard_deviation'] / absolute_sum,
        'new_square_sample_points': point_count,
        'square_samples': samples,
        'smallest_sampled_response': smallest,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scales', nargs='+', type=int, default=[2, 3, 4, 6, 8, 12, 16])
    parser.add_argument('--means', type=Path,
                        default=Path(__file__).with_name('mixed_sieve_response.json'))
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    mean_bytes = args.means.read_bytes()
    existing_means = {r['m']: r['mean_decimal120']
                      for r in json.loads(mean_bytes)['records']}
    records = []
    for m in args.scales:
        record = diagnose(m, existing_means)
        records.append(record)
        print(json.dumps({'m': m,
                          'variance': record['periodic_energy']['centered_variance'],
                          'std_over_abs_sum': record['standard_deviation_over_absolute_sum'],
                          'square_samples': record['new_square_sample_points'],
                          'max_sampled_residual': max(
                              s['upper_weight_square_residual'] for s in record['square_samples'])},
                         sort_keys=True), flush=True)
    result = {
        'scope': 'New finite energy and reciprocal-square diagnostics on arbitrary clocks.',
        'source_means_sha256': hashlib.sha256(mean_bytes).hexdigest(),
        'energy_formula': 'Var(F)=sum_q J_2(q)*(sum_q|d h(d)/d)^2/12.',
        'residual_formula': 'E(y)=(F(y)-y*max(S,0))_+.',
        'sampling_formula': 'SquareCost(x)=sum_p^2<x/U log(p)*E(x/p^2)/x.',
        'boundaries': [
            'No old response grid is rescanned.',
            'No selected N, CA weights delta, or source-phase law is constructed.',
            'Period energy does not certify these sampled costs or a uniform supremum.',
            'Only five x values per scale and the square layer are sampled.',
            'The existing high-precision means are reused; energy is floating point.',
            'Decimal spot agreement is neither interval arithmetic nor Lean verification.',
            'No asymptotic rate, MC_eta, Robin improvement or RH conclusion is established.',
        ],
        'records': records,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
