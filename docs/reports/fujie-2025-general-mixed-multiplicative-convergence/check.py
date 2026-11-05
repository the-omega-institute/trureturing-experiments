#!/usr/bin/env python3
"""Exact finite checks for Fujie's general mixed convergence candidate.

These tests establish finite examples, not an asymptotic convergence theorem.
Run with Python 3 and SymPy. No floating-point roots or tolerances are used.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, permutations
from math import comb, factorial
from random import Random
import sympy as sp

X = sp.Symbol('x')
SEED = 20261105


def elementary(roots):
    e = [Fraction(1)] + [Fraction(0)] * len(roots)
    for r in roots:
        for k in range(len(roots), 0, -1):
            e[k] += Fraction(r) * e[k - 1]
    return e


@lru_cache(maxsize=None)
def convolution(a, b):
    assert len(a) == len(b) and len(a) > 0
    n = len(a)
    ea, eb = elementary(a), elementary(b)
    return sp.Poly.from_list([
        sp.Rational((-1)**k * ea[k] * eb[k] / comb(n, k))
        for k in range(n + 1)], X, domain=sp.QQ)


def root_polynomial(a):
    return sp.Poly.from_list([
        sp.Rational((-1)**k * e) for k, e in enumerate(elementary(a))
    ], X, domain=sp.QQ)


def multiplicity(f, z):
    k = 0
    while f.eval(z) == 0:
        k += 1
        f = f.exquo(sp.Poly(X - z, X, domain=sp.QQ))
    return k


def count_in_interval(f, interval):
    lo, hi = interval
    if lo == hi:
        return multiplicity(f, lo)
    # All factors are square-free; multiplicities are restored explicitly.
    assert f.eval(lo) != 0 and f.eval(hi) != 0
    return sum(m * int(g.count_roots(lo, hi)) for g, m in f.sqf_list()[1])


def exact_cdf_count_difference(f, g):
    n = f.degree()
    assert g.degree() == n and f.LC() == g.LC() == 1
    intervals = (f * g).sqf_part().intervals(eps=sp.Rational(1, 10**6))
    # Sturm isolation must certify all roots real, including multiplicities.
    assert sum(count_in_interval(f, i) for i, _ in intervals) == n
    assert sum(count_in_interval(g, i) for i, _ in intervals) == n
    delta = 0
    worst = 0
    for interval, _ in intervals:
        delta += count_in_interval(f, interval) - count_in_interval(g, interval)
        worst = max(worst, abs(delta))
    assert delta == 0
    return worst


def check_permutation_formula(a, b, f):
    n = len(a)
    total = sp.Poly(0, X, domain=sp.QQ)
    for perm in permutations(b):
        total += root_polynomial(tuple(x*y for x, y in zip(a, perm)))
    assert total.mul_ground(sp.Rational(1, factorial(n))) == f


def main():
    print('Exact arithmetic: Fraction coefficients; SymPy QQ polynomials and Sturm isolation')
    print('SymPy version:', sp.__version__)
    print('Random seed:', SEED)
    fidelity = 0
    small = 0
    for n in range(1, 5):
        for a in combinations_with_replacement((-1, 0, 1), n):
            p = root_polynomial(a)
            # Source Section 3.1: (x(x-1)^(n-1)) boxtimes_n p = x p\'/n.
            qroots = (0,) + (1,) * (n - 1)
            assert convolution(a, qroots) == sp.Poly(X, X) * p.diff() / n
            fidelity += 1
            for b in combinations_with_replacement((0, 1, 2), n):
                f = convolution(a, b)
                check_permutation_formula(a, b, f)
                for alpha in set(a) - {0}:
                    for beta in set(b) - {0}:
                        lower = max(a.count(alpha) + b.count(beta) - n, 0)
                        assert multiplicity(f, alpha*beta) >= lower
                assert multiplicity(f, 0) >= max(a.count(0), b.count(0))
                small += 1
    print('Fidelity: source differential-operator identity, degrees 1..4:', fidelity, 'cases')
    print('Permutation identity and atom lower bounds, all multisets on {-1,0,1} x {0,1,2}, degrees 1..4:', small, 'pairs')
    rng = Random(SEED)
    random_count = 0
    edge_count = 0
    worst = 0
    for n in range(1, 13):
        cases = []
        for _ in range(6):
            a = tuple(Fraction(rng.randint(-9, 9), rng.randint(1, 3)) for _ in range(n))
            b = tuple(Fraction(rng.randint(0, 9), rng.randint(1, 3)) for _ in range(n))
            cases.append((a, b, True))
        cases.extend([
            ((0,) * n, (0,) * n, False),
            ((-1,) * n, (1,) * n, False),
            (tuple((-1)**j for j in range(n)), tuple(j % 2 for j in range(n)), False),
        ])
        for a, b, is_random in cases:
            f = convolution(a, b)
            i = rng.randrange(n)
            ar = list(a)
            ar[i] = Fraction(rng.randint(-12, 12), rng.randint(1, 3))
            br = list(b)
            br[i] = Fraction(rng.randint(0, 12), rng.randint(1, 3))
            for other in (convolution(tuple(ar), b), convolution(a, tuple(br))):
                discrepancy = exact_cdf_count_difference(f, other)
                assert discrepancy <= 1, (n, a, b, discrepancy)
                worst = max(worst, discrepancy)
                random_count += int(is_random)
                edge_count += int(not is_random)
    print('One-root replacement: random exact rational inputs, degrees 1..12:', random_count, 'comparisons')
    print('One-root replacement: zero/repeated/mixed-sign boundary inputs, degrees 1..12:', edge_count, 'comparisons')
    print('Maximum unnormalized CDF discrepancy:', worst, '(bound 1, hence d_K <= 1/n)')
    # This is a warning against using a stronger zero equality; not a refutation of the target.
    assert convolution((-1, 1), (0, 1)) == sp.Poly(X**2, X, domain=sp.QQ)
    print('Zero-atom pitfall: (x^2-1) boxtimes_2 x(x-1) = x^2; only the lower bound is used')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
