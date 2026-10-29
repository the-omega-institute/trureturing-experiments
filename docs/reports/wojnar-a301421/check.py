#!/usr/bin/env python3
"""Exact finite checks for Wojnar's OEIS A301421 GF conjecture.

Sources: arXiv:1706.08381v2, equation (2) and Table GW.n=6;
OEIS A301421 revision 31, NAME, FORMULA and DATA.
Only integers and fractions.Fraction are used. Finite agreement is not a proof.
"""
from fractions import Fraction
from math import comb, factorial, prod

N = 6
MAX_DEGREE = 24
ZERO = (0,) * N
SOURCE_TERMS = [
    1, 6, 46, 371, 3026, 24707, 201748, 1647429, 13452565,
    109850886, 897019828, 7324880157, 59813470848, 488424550081,
    3988374821616, 32568251770049, 265945672309613, 2171657880797162,
    17733313387923690, 144806604435722311, 1182461068019218530,
    9655734852907204771,
]


def partition_vectors(k, i=1):
    if i == N:
        if k % i == 0:
            yield (k // i,)
        return
    for multiplicity in range(k // i + 1):
        for tail in partition_vectors(k - i * multiplicity, i + 1):
            yield (multiplicity,) + tail


def literal_polynomial(k):
    terms = {}
    for m in partition_vectors(k):
        r = sum(m)
        assert r > 0 and sum(i * v for i, v in enumerate(m, 1)) == k
        coefficient = (
            Fraction(k, N * r)
            * (-1) ** (k + r)
            * Fraction(factorial(r), prod(factorial(v) for v in m))
            * prod(comb(N, i) ** v for i, v in enumerate(m, 1))
        )
        assert m not in terms
        terms[m] = coefficient
    return terms


def add_term(poly, m, c):
    value = poly.get(m, Fraction(0)) + c
    if value:
        poly[m] = value
    else:
        poly.pop(m, None)


def newton_polynomials(limit):
    # Independently construct unnormalized power sums from Newton's identity:
    # p_k = -sum_{1 <= i < k} d_i p_{k-i} - k d_k,
    # d_i = (-1)^i binomial(6,i) t_i for 1 <= i <= 6, and zero otherwise.
    power_sums = [{ZERO: Fraction(N)}]
    for k in range(1, limit + 1):
        poly = {}
        for i in range(1, min(N, k - 1) + 1):
            scalar = -((-1) ** i * comb(N, i))
            for m, c in power_sums[k - i].items():
                shifted = list(m)
                shifted[i - 1] += 1
                add_term(poly, tuple(shifted), scalar * c)
        if k <= N:
            m = list(ZERO)
            m[k - 1] = 1
            add_term(poly, tuple(m), Fraction(-k * (-1) ** k * comb(N, k)))
        power_sums.append(poly)
    return [{m: c / N for m, c in p.items()} for p in power_sums]


def gf_coefficients(limit):
    # Coefficients of the printed numerator and denominator, constant term first.
    numerator = [1, -1, -5, -10, -10, -5, -1]
    denominator = [1, -7, -9, -5, 5, 9, 5, 1]
    # Check both factorizations by independent integer polynomial multiplication.
    b = [comb(6, i) for i in range(7)]
    q = [2 - b[0]] + [-x for x in b[1:]]
    expanded = [0] * 8
    for i, c in enumerate(q):
        expanded[i] += c
        expanded[i + 1] -= c
    assert expanded == denominator and denominator[0] == 1
    assert numerator == [1] + [-comb(5, i) for i in range(6)]
    values = []
    for n in range(limit):
        rhs = Fraction(numerator[n] if n < len(numerator) else 0)
        values.append(rhs - sum(
            denominator[i] * values[n - i]
            for i in range(1, min(n, len(denominator) - 1) + 1)
        ))
    return values


def elementary_means(data):
    elementary = [Fraction(1)] + [Fraction(0)] * N
    for x in data:
        for j in range(N, 0, -1):
            elementary[j] += x * elementary[j - 1]
    return [elementary[i] / comb(N, i) for i in range(1, N + 1)]


def evaluate(poly, t):
    return sum((c * prod(x ** v for x, v in zip(t, m))
                for m, c in poly.items()), Fraction(0))


def main():
    literal = [literal_polynomial(k) for k in range(1, MAX_DEGREE + 1)]
    newton = newton_polynomials(MAX_DEGREE)
    assert all(literal[k - 1] == newton[k] for k in range(1, MAX_DEGREE + 1))
    # Source Table GW.n=6 gives these actual polynomials, not just their masses.
    assert literal[0] == {(1, 0, 0, 0, 0, 0): Fraction(1)}
    assert literal[1] == {(2, 0, 0, 0, 0, 0): Fraction(6),
                          (0, 1, 0, 0, 0, 0): Fraction(-5)}
    assert literal[2] == {(3, 0, 0, 0, 0, 0): Fraction(36),
                          (1, 1, 0, 0, 0, 0): Fraction(-45),
                          (0, 0, 1, 0, 0, 0): Fraction(10)}
    masses = [sum((c for c in p.values() if c > 0), Fraction(0)) for p in literal]
    assert masses[:len(SOURCE_TERMS)] == SOURCE_TERMS
    assert masses[:7] == [1, 6, 46, 371, 3026, 24707, 201748]
    assert masses == gf_coefficients(MAX_DEGREE)
    # Exact evaluations include repeated roots, zero data, negatives and fractions.
    data_sets = [
        tuple(map(Fraction, (0, 0, 0, 0, 0, 0))),
        tuple(map(Fraction, (1, 1, 1, 1, 1, 1))),
        tuple(map(Fraction, (-3, -1, 0, 1, 2, 5))),
        (Fraction(-1, 2), Fraction(2, 3), Fraction(3, 5),
         Fraction(-7, 4), Fraction(0), Fraction(2)),
    ]
    for data in data_sets:
        t = elementary_means(data)
        for k, poly in enumerate(literal, 1):
            assert evaluate(poly, t) == sum(x ** k for x in data) / N
    for mass, poly in zip(masses, literal):
        assert sum(poly.values()) == 1
        assert 2 * mass == sum(abs(c) for c in poly.values()) + 1
    print('FIDELITY: exact Table GW.n=6 polynomials at k=1,2,3 and positive masses at k=1..7 reproduced')
    print('OEIS: all 22 DATA terms in A301421 revision 31 reproduced')
    print('INDEPENDENCE: partition coefficients equal aggregated Newton-recursion polynomials for k=1..24')
    print('SOURCE_BRIDGE: 4 exact six-data sets, k=1..24; repeated, zero, negative and rational data included')
    print('GF: exact coefficients agree for k=1..24; denominator constant term = 1')
    print('MASS: signed mass = 1 and twice positive mass = unsigned mass + 1 for k=1..24')
    print('VALUES k=1..24: ' + ','.join(str(x) for x in masses))
    print('MONOMIAL_COUNTS k=1..24: ' + ','.join(str(len(p)) for p in literal))
    print('SCOPE: finite checks only; no all-degree proof and no Lean proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
