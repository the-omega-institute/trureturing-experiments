#!/usr/bin/env python3
"""Exact finite coefficient and lattice checks; no universal proof is certified.

Source: Daniels--Huber--McLaughlin--Ye, arXiv:2606.06701v1,
Theorem 1.2, Example 4.1(1), and the first experimental example in Section 5.
Companion: arXiv:2603.04666v1, Section 6.
SPDX-License-Identifier: Apache-2.0
"""

from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import gcd, isqrt


def q_product_factors(a, m, cutoff):
    """All factors of Q(q^a,q^m) that can affect degrees <= cutoff."""
    assert 0 < 2 * a < m
    return [
        e
        for start, step in ((a, m), (m - a, m), (m, m),
                            (m + 2 * a, 2 * m), (m - 2 * a, 2 * m))
        for e in range(start, cutoff + 1, step)
    ]


def product_coefficients(exponents, cutoff):
    c = [1] + [0] * cutoff
    for e in exponents:
        assert e > 0
        for d in range(cutoff, e - 1, -1):
            c[d] -= c[d - e]
    return c


def quintuple_terms(a, m, cutoff):
    """Signed terms of the source's bilateral identity, with their lattice x.

    Completing the square bounds |6mn+6a+-m| by the integer square root.
    The generous n range below follows directly from this bound, so no term
    affecting a checked coefficient is omitted.
    """
    radius = isqrt(24 * m * cutoff + (6 * a - m) ** 2)
    bound = (radius + 6 * a + m) // (6 * m) + 2
    terms = []
    for n in range(-bound, bound + 1):
        branches = (
            (m * (3 * n * n - n) // 2 + 3 * a * n,
             1, 6 * m * n + 6 * a - m),
            (m * (3 * n * n + n) // 2 + a * (3 * n + 1),
             -1, 6 * m * n + 6 * a + m),
        )
        for exponent, sign, x in branches:
            assert exponent >= 0
            assert x * x - (6 * a - m) ** 2 == 24 * m * exponent
            if exponent <= cutoff:
                terms.append((exponent, sign, x))
    return terms


def theta_coefficients(parameters, m, cutoff):
    c = {0: 1}
    for a in parameters:
        out = defaultdict(int)
        for exponent, sign, _ in quintuple_terms(a, m, cutoff):
            for d, value in c.items():
                if d + exponent <= cutoff:
                    out[d + exponent] += value * sign
        c = out
    return [c.get(d, 0) for d in range(cutoff + 1)]


def checked_product(parameters, m, cutoff):
    factors = [e for a in parameters for e in q_product_factors(a, m, cutoff)]
    direct = product_coefficients(factors, cutoff)
    independent = theta_coefficients(parameters, m, cutoff)
    assert direct == independent, (parameters, m)
    return direct


def check_progressions(c, m, residues, last_t):
    for t in range(last_t + 1):
        for residue in residues:
            assert c[m * t + residue] == 0, (m, t, residue)


def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3))
             for j in range(3)] for i in range(3)]


def transpose(a):
    return [list(x) for x in zip(*a)]


def lattice_checks(cutoff):
    r = (1, 5, 6)
    ident = [[int(i == j) for j in range(3)] for i in range(3)]
    reflection = [[Fraction(ident[i][j]) - Fraction(r[i] * r[j], 31)
                   for j in range(3)] for i in range(3)]
    cyclic_matrix = [[0, 1, 0], [0, 0, -1], [-1, 0, 0]]
    assert matmul(reflection, reflection) == ident
    assert matmul(transpose(reflection), reflection) == ident
    assert matmul(transpose(cyclic_matrix), cyclic_matrix) == ident
    assert matmul(matmul(cyclic_matrix, cyclic_matrix), cyclic_matrix) == ident
    for j in range(3):
        p = ident
        for _ in range((-j) % 3):
            p = matmul(cyclic_matrix, p)
        signed = [[-v for v in row] for row in p]
        combined = matmul(signed, reflection)
        assert matmul(combined, combined) == ident

    def cyclic(x):
        return (x[1], -x[2], -x[0])

    def weight(x):
        eps = [1 if v % 6 == 1 else -1 for v in x]
        assert all(v % 6 in (1, 5) for v in x)
        return -eps[0] * eps[1] * eps[2]

    def transform(x, j):
        dot = sum(a * b for a, b in zip(r, x))
        assert dot % 31 == 0
        k = dot // 31
        assert k % 2 == 0
        y = tuple(v - k * a for v, a in zip(x, r))
        assert weight(y) == weight(x)
        for _ in range((-j) % 3):
            y = cyclic(y)
        return tuple(-v for v in y)

    shells = defaultdict(set)
    terms = [quintuple_terms(a, 31, cutoff) for a in r]
    for triple in product(*terms):
        degree = sum(z[0] for z in triple)
        if degree <= cutoff and degree % 31 in (6, 16, 18):
            x = tuple(z[2] for z in triple)
            assert sum(v * v for v in x) == 744 * degree + 651
            assert weight(x) == triple[0][1] * triple[1][1] * triple[2][1]
            shells[degree].add(x)
    point_count = 0
    for degree, shell in sorted(shells.items()):
        j = {6: 0, 18: 1, 16: 2}[degree % 31]
        images = set()
        for x in shell:
            y = transform(x, j)
            assert y in shell
            assert sum(v * v for v in y) == sum(v * v for v in x)
            assert all((v - 6 * a) % 31 == 0 for v, a in zip(y, r))
            assert weight(y) == -weight(x)
            images.add(y)
            assert transform(y, j) == x
        assert images == shell
        assert sum(weight(x) for x in shell) == 0
        point_count += len(shell)
    print('LATTICE finite nonempty shells:', len(shells), 'points:', point_count)
    print('LATTICE all three combined maps square to I over the rationals; signs reverse')


def main():
    # Fidelity: an established two-factor result, separate from the target.
    c = checked_product((6, 1), 37, 37 * 20 + 22)
    check_progressions(c, 37, (19, 22), 20)
    assert c[:8] == [1, -1, 0, 0, 0, 0, -1, 1]
    print('FIDELITY Example 4.1(1), b=1: 42 zeros through t=20; product=quintuple sum')

    cutoff = 31 * 20 + 18
    c = checked_product((1, 5, 6), 31, cutoff)
    check_progressions(c, 31, (6, 16, 18), 20)
    print('TARGET coefficients a_0..a_30:', c[:31])
    print('TARGET 63 zeros, t=0..20, degrees 6..638; all 639 coefficients agree independently')

    for parameters, m, residues in [((1, 3, 4), 13, (2, 4, 10)),
                                     ((2, 3, 5), 19, (4, 5, 16)),
                                     ((3, 5, 8), 49, (8, 10, 46))]:
        d = m * 20 + max(residues)
        v = checked_product(parameters, m, d)
        check_progressions(v, m, residues, 20)
        print('RELATED', parameters, 'modulus', m, 'residues', residues, '63 zeros')

    count = 0
    for a in range(1, 9):
        for b in range(a + 1, 9):
            if gcd(a, b) != 1 or (a - b) % 3 == 0:
                continue
            m = a * a + a * b + b * b
            zeta = b * pow(a, -1, m) % m
            assert (zeta * zeta + zeta + 1) % m == 0
            assert m % 6 == 1 and 2 * (a + b) < m
            residues = tuple(sorted({(a + b + 3 * (pow(zeta, j, m) - 1)) % m
                                     for j in range(3)}))
            d = m * 5 + max(residues)
            v = checked_product((a, b, a + b), m, d)
            check_progressions(v, m, residues, 5)
            count += 1
    print('FAMILY', count, 'coprime pairs 1<=A<B<=8 with A!=B mod 3; t=0..5')
    lattice_checks(cutoff)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
