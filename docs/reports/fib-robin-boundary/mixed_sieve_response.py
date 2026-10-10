#!/usr/bin/env python3
"""Finite diagnostics for the proposed mixed-sieve negative-part estimate.

These are arbitrary real-clock samples A=m**8, not selected Robin sources.
Decimal agreement is a numerical check, not an interval or Lean certificate.
Only the Python standard library is required.
"""

import argparse
from array import array
from collections import defaultdict
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path


def mobius_table(limit):
    mu = [1] * (limit + 1)
    mu[0] = 0
    composite = bytearray(limit + 1)
    for p in range(2, limit + 1):
        if composite[p]:
            continue
        for n in range(p, limit + 1, p):
            composite[n] = 1
            mu[n] = -mu[n]
        for n in range(p * p, limit + 1, p * p):
            mu[n] = 0
    return mu


def float_weights(m, mu):
    log_m = math.log(m)
    lam = [(d, mu[d] * math.log(m / d) / log_m)
           for d in range(1, m) if mu[d]]
    theta = [(e, float(mu[e]) if e <= m * m else
              mu[e] * math.log(m ** 3 / e) / log_m)
             for e in range(1, m ** 3) if mu[e]]
    return lam, theta


def lcm_coefficients(lam, theta):
    terms = defaultdict(list)
    for d, a in lam:
        for e, b in theta:
            terms[d * e // math.gcd(d, e)].append(a * b)
    return sorted((d, math.fsum(values)) for d, values in terms.items())


def grid_response(lam, theta, limit):
    left = array('d', [0.0]) * (limit + 1)
    right = array('d', [0.0]) * (limit + 1)
    for weights, values in [(lam, left), (theta, right)]:
        for d, a in weights:
            for n in range(d, limit + 1, d):
                values[n] += a
    response = array('d', [0.0]) * (limit + 1)
    total = correction = 0.0
    for n in range(1, limit + 1):
        term = left[n] * right[n] - correction
        new_total = total + term
        correction = (new_total - total) - term
        total = new_total
        response[n] = total
    return response


def decimal_response(m, mu, ys, precision):
    """Evaluate the pair sum directly, without grouped LCM coefficients."""
    with localcontext() as ctx:
        ctx.prec = precision
        log_m = Decimal(m).ln()
        lam = [(d, Decimal(mu[d]) * (log_m - Decimal(d).ln()) / log_m)
               for d in range(1, m) if mu[d]]
        theta = [(e, Decimal(mu[e]) if e <= m * m else
                  Decimal(mu[e]) * (3 * log_m - Decimal(e).ln()) / log_m)
                 for e in range(1, m ** 3) if mu[e]]
        mean = Decimal(0)
        responses = {y: Decimal(0) for y in ys}
        for d, a in lam:
            for e, b in theta:
                k = d * e // math.gcd(d, e)
                weight = a * b
                mean += weight / k
                for y in ys:
                    responses[y] += weight * (y // k)
        return mean, responses


def diagnose(m):
    if m < 2:
        raise ValueError('Each scale m must be at least two.')
    if m ** 4 > 1_000_000:
        raise ValueError('The enumerated grid must have at most 1,000,000 points.')
    mu = mobius_table(m ** 3)
    lam, theta = float_weights(m, mu)
    h = lcm_coefficients(lam, theta)
    grid = grid_response(lam, theta, m ** 4)
    minimum_at = min(range(1, len(grid)), key=grid.__getitem__)
    ys = sorted({1, m ** 2, m ** 3, m ** 4, m ** 6, m ** 8, minimum_at})
    mean80, responses80 = decimal_response(m, mu, ys, 80)
    mean120, responses120 = decimal_response(m, mu, ys, 120)
    with localcontext() as ctx:
        ctx.prec = 120
        differences = [abs(mean80 - mean120)] + [
            abs(responses80[y] - responses120[y]) for y in ys]
        precision_difference = max(differences)
        if precision_difference > Decimal('1e-55'):
            raise ArithmeticError('The two Decimal precisions disagree.')
        corrections120 = {y: responses120[y] - y * mean120 for y in ys}
    floor_values = {y: math.fsum(v * (y // d) for d, v in h) for y in ys}
    mismatch = max(abs(grid[y] - floor_values[y]) for y in ys if y < len(grid))
    if mismatch > 1e-9:
        raise ArithmeticError('Direct divisor and grouped floor evaluations disagree.')
    numerical_tolerance = 1e-10
    observed_negative = sum(grid[y] < -numerical_tolerance for y in range(1, len(grid)))
    initial_error = max(abs(grid[y] - 1) for y in range(1, m ** 2 + 1))
    if initial_error > numerical_tolerance:
        raise ArithmeticError('The response below the first cutoff disagrees with one.')
    mean_float = math.fsum(v / d for d, v in h)
    return {
        'm': m, 'A': m ** 8, 'R': m, 'U': m ** 2, 'U1': m ** 3,
        'support_upper_bound': m ** 4,
        'enumerated_grid': [1, m ** 4],
        'grouped_lcm_count': len(h),
        'mean_float': mean_float,
        'mean_decimal120': str(mean120),
        'sum_absolute_grouped_coefficients_float': math.fsum(abs(v) for _, v in h),
        'observed_negative_grid_points_below_minus_1e_minus_10': observed_negative,
        'grid_minimum': {'y': minimum_at, 'float': grid[minimum_at],
                         'decimal120': str(responses120[minimum_at])},
        'initial_cutoff_maximum_numerical_error': initial_error,
        'selected_direct_vs_grouped_maximum_error': mismatch,
        'decimal80_vs_decimal120_maximum_error': str(precision_difference),
        'selected_responses': [
            {'y': y, 'H_decimal120': str(responses120[y]),
             'floor_correction_decimal120': str(corrections120[y])}
            for y in ys],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scales', nargs='+', type=int, default=[2, 3, 4, 6, 8, 12, 16])
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    records = []
    for m in args.scales:
        record = diagnose(m)
        records.append(record)
        print(json.dumps({k: record[k] for k in [
            'm', 'mean_float', 'grid_minimum',
            'observed_negative_grid_points_below_minus_1e_minus_10']},
            sort_keys=True), flush=True)
    result = {
        'scope': 'Finite coefficient diagnostics, not actual selected integer sources.',
        'weights': 'rho_z(d)=mu(d)log(z/d)1[d<=z]; lambda=rho_R/log(R); '
                   'theta_prime=(rho_U1-rho_U)/log(U1/U).',
        'response': 'H(y)=sum_n<=y (sum_d|n lambda(d))(sum_e|n theta_prime(e)).',
        'mean': 'S=sum_[d,e] lambda(d)theta_prime(e)/lcm(d,e).',
        'boundaries': [
            'No actual N or selected CA/GA1/right-tail-max source is constructed.',
            'MC_eta is not measured or proved.',
            'Neither grid signs nor Decimal agreement are a uniform asymptotic estimate.',
            'The computed corrections retain floor factors at and across every cutoff.',
            'No interval arithmetic or Lean kernel verification is performed.',
        ],
        'records': records,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
