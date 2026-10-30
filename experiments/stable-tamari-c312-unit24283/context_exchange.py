#!/usr/bin/env python3
"""Exact bounded test of the Thm 4.17 two-summand exchange.
Standard library only. No enumeration of top vectors b. Stops at first failure.
Exponent tuples list x_1,...,x_n,t; a denominator factor is 1 - monomial.
Source formula:  arXiv:2610.39257v1,
The Stable Tamari Lattice (Chau et al.), CC BY 4.0. This file is original code.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent


def valid(h):
    if any(not 1 <= left <= right for right, left in enumerate(h, 1)):
        return False
    intervals = [(left, right) for right, left in enumerate(h, 1)]
    for i, (a, b) in enumerate(intervals):
        for c, d in intervals[i + 1:]:
            disjoint = b < c or d < a
            contains = (a <= c and d <= b) or (c <= a and b <= d)
            if not (disjoint or contains):
                return False
    return True


def eligible(h0, p):
    n = len(h0)
    if not (h0[p - 1] == p and h0[p] == p + 1):
        return None
    h1 = h0[:p] + (p,) + h0[p + 1:]
    if not valid(h1):
        return None
    for r, left in enumerate(h0, 1):
        if r in (p, p + 1):
            continue
        disjoint = r < p or left > p + 1
        contains = left <= p and r >= p + 1
        if not (disjoint or contains):
            return None
    return h1


def interval_exp(n, left, right, colored=False):
    return tuple(int(left <= i <= right) for i in range(1, n + 1)) + (int(colored),)


def summand(h):
    n = len(h)
    a = tuple(sum(left <= i <= r for r, left in enumerate(h, 1)) - 1
              for i in range(1, n + 1))
    assert min(a) >= 0
    uncolored = []
    u_data = []
    for j in range(1, n + 1):
        if h[j - 1] != j:
            continue
        previous = a[j - 2] if j > 1 else 0
        for gamma in range(previous, a[j - 1] + 1):
            k = j
            while k < n and a[k] >= gamma:
                k += 1
            assert all(a[q - 1] >= gamma for q in range(j, k + 1))
            assert k == n or a[k] < gamma
            uncolored.append(interval_exp(n, j, k))
            u_data.append({'j': j, 'gamma': gamma, 'k': k})
    colored = [interval_exp(n, left, r, True) for r, left in enumerate(h, 1)]
    return {'h': h, 'a': a, 'numerator': a + (0,),
            'U': u_data, 'uncolored': uncolored, 'colored': colored,
            'denominator': uncolored + colored}


def swap(exp, p):
    result = list(exp)
    result[p - 1], result[p] = result[p], result[p - 1]
    return tuple(result)


def add_term(poly, exp, coefficient):
    poly[exp] = poly.get(exp, 0) + coefficient
    if poly[exp] == 0:
        del poly[exp]


def binomial_multiply(poly, factor):
    result = dict(poly)
    for exp, coefficient in poly.items():
        shifted = tuple(a + b for a, b in zip(exp, factor))
        add_term(result, shifted, -coefficient)
    return result


def clear_difference(s0, s1, p):
    terms = []
    for s in (s0, s1):
        terms.append((1, s['numerator'], Counter(s['denominator'])))
    for s in (s0, s1):
        terms.append((-1, swap(s['numerator'], p),
                      Counter(swap(e, p) for e in s['denominator'])))
    common = Counter()
    for _, _, factors in terms:
        common |= factors
    result = {}
    for sign, numerator, factors in terms:
        poly = {numerator: sign}
        for factor, multiplicity in (common - factors).items():
            for _ in range(multiplicity):
                poly = binomial_multiply(poly, factor)
        for exp, coefficient in poly.items():
            add_term(result, exp, coefficient)
    return result, common


def monomial_value(exp, values):
    result = Fraction(1)
    for e, value in zip(exp, values):
        result *= value ** e
    return result


def rational_value(s, values, p=None):
    transform = (lambda e: swap(e, p)) if p is not None else (lambda e: e)
    numerator = monomial_value(transform(s['numerator']), values)
    denominator = Fraction(1)
    for e in s['denominator']:
        factor = 1 - monomial_value(transform(e), values)
        assert factor != 0
        denominator *= factor
    return numerator / denominator


def polynomial_value(poly, values):
    return sum((coefficient * monomial_value(exp, values)
                for exp, coefficient in poly.items()), Fraction(0))


def main():
    inspections = []
    witness = None
    # Bounded 2 <= n <= 4 domain, with first-failure exit.
    for n in range(2, 5):
        for h0 in product(*(range(1, r + 1) for r in range(1, n + 1))):
            if not valid(h0):
                continue
            for p in range(1, n):
                h1 = eligible(h0, p)
                if h1 is None:
                    continue
                s0, s1 = summand(h0), summand(h1)
                numerator, common = clear_difference(s0, s1, p)
                inspections.append({'n': n, 'p': p, 'h0': h0, 'h1': h1,
                                    'cleared_difference_zero': not numerator})
                if numerator:
                    witness = (n, p, h0, h1, s0, s1, numerator, common)
                    break
            if witness:
                break
        if witness:
            break
    assert witness is not None, 'No obstruction: report missing generic bridge.'
    n, p, h0, h1, s0, s1, numerator, common = witness
    assert (n, p, h0, h1) == (3, 1, (1, 2, 3), (1, 1, 3))
    assert s0['numerator'] == (0, 0, 0, 0)
    assert s1['numerator'] == (1, 0, 0, 0)
    assert Counter(s0['denominator']) == Counter([
        (1,1,1,0), (0,1,1,0), (0,0,1,0),
        (1,0,0,1), (0,1,0,1), (0,0,1,1)])
    assert Counter(s1['denominator']) == Counter([
        (1,1,1,0), (1,0,0,0), (0,0,1,0),
        (1,0,0,1), (1,1,0,1), (0,0,1,1)])
    # Verify full symbolic certificate Q*D=(x-y)(1-xyz)(1-z)(1-txy).
    expected = {(1,0,0,0): 1, (0,1,0,0): -1}
    canceled = [(1,1,1,0), (0,0,1,0), (1,1,0,1)]
    for factor in canceled:
        expected = binomial_multiply(expected, factor)
    assert numerator == expected
    remaining = common - Counter(canceled)
    assert remaining == Counter([
        (1,0,0,1), (0,1,0,1), (0,0,1,1),
        (1,0,0,0), (0,1,0,0), (1,0,1,0), (0,1,1,0)])
    values = (Fraction(1,5), Fraction(1,7), Fraction(1,11), Fraction(2))
    lhs = rational_value(s0, values) + rational_value(s1, values)
    rhs = rational_value(s0, values, p) + rational_value(s1, values, p)
    difference = lhs - rhs
    factored_difference = (values[0] - values[1])
    for factor, multiplicity in remaining.items():
        factored_difference /= (1 - monomial_value(factor, values)) ** multiplicity
    assert difference == factored_difference and difference != 0
    common_value = Fraction(1)
    for factor, multiplicity in common.items():
        common_value *= (1 - monomial_value(factor, values)) ** multiplicity
    assert polynomial_value(numerator, values) == common_value * difference
    data = {
        'purpose': 'Exact finite obstruction to two-summand context exchange; not a counterexample to Conjecture 3.12',
        'source': {'title': 'The Stable Tamari Lattice', 'arxiv': '2610.39257v1',
                   'license': 'CC BY 4.0', 'formula': 'Theorem 4.17'},
        'domain': {'n_min': 2, 'n_max': 4, 'stop': 'first exact nonzero difference',
                   'order': 'n, lexicographic h0, ascending p', 'top_b_search': False},
        'inspections': inspections, 'n4_inspected': False,
        'witness': {'n': n, 'p': p, 'h0': h0, 'h1': h1,
                    's0': s0, 's1': s1,
                    'D_definition': 'W_h0(x;t)+W_h1(x;t)-W_h0(sigma_p x;t)-W_h1(sigma_p x;t)',
                    'common_denominator_factors': [
                        {'exponents': e, 'multiplicity': m} for e,m in sorted(common.items())],
                    'cleared_numerator': [
                        {'exponents': e, 'coefficient': c} for e,c in sorted(numerator.items())],
                    'cleared_numerator_factorization': '(x-y)(1-xyz)(1-z)(1-txy)',
                    'reduced_difference': '(x-y)/((1-tx)(1-ty)(1-tz)(1-x)(1-y)(1-xz)(1-yz))',
                    'numeric_check': {'x': '1/5', 'y': '1/7', 'z': '1/11', 't': '2',
                                      'lhs': str(lhs), 'rhs': str(rhs),
                                      'lhs_minus_rhs': str(difference)}},
        'verification': {'arithmetic': 'Python integers and fractions.Fraction',
                         'full_sparse_polynomial_identity': True,
                         'Lean_kernel': False},
    }
    path = ROOT / 'context_exchange_data.json'
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'inspections': inspections, 'witness': data['witness']['reduced_difference'],
                      'numeric_difference': str(difference)}, indent=2))


if __name__ == '__main__':
    main()
