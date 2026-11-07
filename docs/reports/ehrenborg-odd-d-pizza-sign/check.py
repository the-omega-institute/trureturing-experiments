#!/usr/bin/env python3
"""Exact finite checks for Ehrenborg, arXiv:2410.13593v1, Conjecture 5.2.

Only integers, Fraction, and exact SymPy polynomials are used as evidence.
This program is a normalization and finite-case check, not a universal proof.
"""
from fractions import Fraction as Q
from functools import cache
from itertools import combinations
from math import comb, factorial
import sympy as sp


@cache
def h(r):
    """Gamma(1/2, scale 4) integer moment."""
    return factorial(2 * r) // factorial(r)


@cache
def direct_pair(p, q):
    """[a^(2p)b^(2q)z^(p+q-1)] E_z((a+b)/sqrt2)E_z((a-b)/sqrt2).

    Expand both odd powers directly; sqrt2 cancels to an integer power of 2.
    No moment formula or proposed coefficient recurrence is used here.
    """
    if p + q == 0:
        return Q(0)
    total = Q(0)
    for r in range(p + q):
        s = p + q - 1 - r
        u, v = 2 * r + 1, 2 * s + 1
        b_coefficient = sum(
            comb(u, j) * comb(v, 2 * q - j) * (-1) ** (2 * q - j)
            for j in range(u + 1) if 0 <= 2 * q - j <= v
        )
        total += Q(b_coefficient,
                   2 ** (p + q) * factorial(r) * factorial(s) * u * v)
    return total


@cache
def moment_pair(p, q):
    """Integral of (x^p y^q - x^q y^p)/(x-y), after polynomial cancellation."""
    if p == q:
        return 0
    if p < q:
        return -moment_pair(q, p)
    return sum(h(r) * h(p + q - 1 - r) for r in range(q, p))


def pfaffian(a):
    """Signed matching expansion, pairing the first remaining index."""
    n = len(a)
    assert n % 2 == 0
    if n == 0:
        return 1
    return sum((-1) ** (j + 1) * a[0][j] * pfaffian([
        [a[r][s] for s in range(1, n) if s != j]
        for r in range(1, n) if r != j
    ]) for j in range(1, n))


def c(m):
    value = Q((-1) ** m, factorial(m))
    for j in range(m):
        value *= Q(1, 2) - j
    return value


def partition_vectors(total, length, upper=None):
    if length == 0:
        if total == 0:
            yield ()
        return
    if upper is None:
        upper = total
    for a in range(min(upper, total), -1, -1):
        for tail in partition_vectors(total - a, length - 1, a):
            yield (a,) + tail


def source_coefficient(n, lam, pair=direct_pair):
    """T_d's multinomial sum factors by disjoint pairs before signed summation."""
    size = n - 1
    assert n >= 3 and n % 2 == 1 and len(lam) == size
    assert all(lam[i] >= lam[i + 1] for i in range(size - 1))
    p = [lam[i] + size - i for i in range(size)]
    d = n * size + 2 * sum(lam)
    m = (d - size) // 2
    a = [[pair(p[i], p[j]) for j in range(size)] for i in range(size)]
    return c(m) * factorial(m) * pfaffian(a)


def moment_coefficient(n, lam):
    size = n - 1
    p = [lam[i] + size - i for i in range(size)]
    normalized = pfaffian([[moment_pair(i, j) for j in p] for i in p])
    m = (n * size + 2 * sum(lam) - size) // 2
    return c(m) * factorial(m) * Q(normalized,
        sp.prod(factorial(2 * i) for i in p))


def exact_moment(poly, xs):
    """Integrate a polynomial against the product moment functional, exactly."""
    p = sp.Poly(poly, *xs)
    return sum(coef * sp.prod(h(e) for e in powers)
               for powers, coef in p.terms())


def main():
    pair_count = 0
    for p in range(13):
        for q in range(13):
            actual = factorial(2 * p) * factorial(2 * q) * direct_pair(p, q)
            assert actual == moment_pair(p, q), (p, q, actual)
            pair_count += 1
    print(f'PAIR_IDENTITY_OK p,q=0..12 cases={pair_count}')

    # Table 2, both columns, for D3, D5, D7. The leading monomial of J
    # has coefficient 2^(-n(n-1)/2); the first Schur polynomial is sum a_i^2.
    expected = {
        3: (Q(-1, 6), Q(-1, 10)),
        5: (Q(-143, 40), Q(-143, 12)),
        7: (Q(-955049953, 336), Q(-5 * 11 * 17 * 19 * 23 * 29 * 31, 16)),
    }
    for n, targets in expected.items():
        for extra, target in enumerate(targets):
            lam = (extra,) + (0,) * (n - 2)
            actual = 2 ** (n * (n - 1) // 2) * source_coefficient(n, lam)
            assert actual == target, (n, extra, actual, target)
            print(f'FIDELITY D{n} degree={n*(n-1)+2*extra} Z_d/J={actual}')

    count = 0
    for n, max_weight in [(3, 12), (5, 6), (7, 3), (9, 2)]:
        local_count = 0
        for total in range(max_weight + 1):
            for lam in partition_vectors(total, n - 1):
                direct = source_coefficient(n, lam)
                moment = moment_coefficient(n, lam)
                assert direct == moment and direct < 0, (n, lam, direct, moment)
                local_count += 1
        print(f'CLAIM_SMALL_CASES_OK n={n} partition_weight=0..{max_weight} cases={local_count}')
        count += local_count
    print(f'COEFFICIENT_CASES_TOTAL {count}')

    # Proposition 5.3 needs the projection onto monomials containing both variables.
    a, b = sp.symbols('a b')
    def td6(u, v):
        return -sp.Rational(1, 8) * (u**5*v/5 + 2*u**3*v**3/9 + u*v**5/5)
    y_poly = sp.expand(td6((a+b)/sp.sqrt(2), (a-b)/sp.sqrt(2)))
    z_special = sp.expand(y_poly - td6(a/sp.sqrt(2), a/sp.sqrt(2))
                         + td6(b/sp.sqrt(2), b/sp.sqrt(2)))
    assert sp.expand(z_special-y_poly) == sp.Rational(7, 720)*(a**6-b**6)
    def all_variables(poly):
        return sum(coef*a**i*b**j for (i, j), coef in sp.Poly(poly, a, b).terms()
                   if i > 0 and j > 0)
    assert all_variables(z_special) == all_variables(y_poly)
    print('PROJECTION_OK n=3 d=6 Z(a,b,0)-Y=7*(a^6-b^6)/720')

    sign_count = 0
    for size in [2, 4, 6]:
        for xs in combinations(range(1, 9), size):
            cauchy = [[Q(0) if i == j else Q(1, xs[i] - xs[j])
                       for j in range(size)] for i in range(size)]
            assert (-1) ** (size // 2) * pfaffian(cauchy) > 0
            sign_count += 1
    print(f'RECIPROCAL_PFAFFIAN_SIGN_OK subsets_of_1..8 sizes=2,4,6 cases={sign_count}')

    # Cancel the full determinant times Pfaffian before applying any moments.
    # This never integrates an individual reciprocal-distance singular term.
    for exponents in [(2, 1), (4, 3, 2, 1), (5, 3, 1, 0), (6, 4, 2, 0)]:
        size = len(exponents)
        xs = sp.symbols(f'x0:{size}')
        det = sp.Matrix([[x ** p for x in xs] for p in exponents]).det()
        kernel = [[0 if i == j else 1 / (xs[i] - xs[j])
                   for j in range(size)] for i in range(size)]
        polynomial = sp.cancel(det * pfaffian(kernel))
        assert sp.denom(polynomial) == 1
        lhs = exact_moment(polynomial, xs)
        rhs = Q(factorial(size), 2 ** (size // 2)) * pfaffian([
            [moment_pair(p, q) for q in exponents] for p in exponents])
        assert lhs == rhs and lhs > 0, (exponents, lhs, rhs)
        print(f'CANCELLATION_MOMENT_IDENTITY_OK exponents={exponents} value={lhs}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
