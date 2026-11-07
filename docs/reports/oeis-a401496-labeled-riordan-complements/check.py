#!/usr/bin/env python3
"""Exact finite diagnostics for OEIS A401496; no claim of an all-order proof.

GF(2)[t]/(t^N) polynomials are integers with coefficient k at bit k.
Enumerate every (g,f) with f(0)=0, build the source adjacency matrix,
and deduplicate graphs before testing their literal fixed-label complements.
The proposed parameter classification is constructed independently and compared
as a set of graphs, including all its nonunique presentations.
"""
from collections import defaultdict
from fractions import Fraction
from math import comb

OEIS = (1, 2, 4, 12, 28, 92, 220, 732, 1756, 5852, 14044, 46812,
        112348, 374492)
A113836 = (1, 3, 11, 27, 91, 219, 731, 1755, 5851, 14043, 46811,
           112347, 374491)
MAX_ENUM = 10


def mul(a, b, N):
    """Carryless product, reduced modulo t^N."""
    c = 0
    while b:
        if b & 1:
            c ^= a
        a <<= 1
        b >>= 1
    return c & ((1 << N) - 1)


def inverse_unit(a, N):
    assert a & 1
    # Independent coefficient recursion for a*b = 1.
    b = 1
    for k in range(1, N):
        parity = sum(((a >> i) & 1) * ((b >> (k - i)) & 1)
                     for i in range(1, k + 1)) % 2
        b |= parity << k
    assert mul(a, b, N) == 1
    return b


def graph(n, g, f):
    N = n - 1
    assert f & 1 == 0
    cols = [g]
    for _ in range(1, N):
        cols.append(mul(cols[-1], f, N))
    mask = 0
    bit = 0
    for i in range(1, n):
        for j in range(i):
            # Source labels I=i+1, J=j+1; [t^(I-2)] g*f^(J-1).
            mask |= ((cols[j] >> (i - 1)) & 1) << bit
            bit += 1
    assert bit == n * (n - 1) // 2
    return mask


def exhaustive(n):
    if n == 1:
        return {0}
    N = n - 1
    return {graph(n, g, f) for g in range(1 << N)
            for f in range(0, 1 << N, 2)}


def classified_units(n):
    N = n - 1
    H = (1 << N) - 1
    t = 2 & H
    result = set()
    per_r = {}
    raw = 0
    for r in range(0, 1 << N, 2):
        g = H ^ r
        invg = inverse_unit(g, N)
        alpha = mul(r, invg, N)
        fixed = set()
        for delta in range(0, 1 << N, 2):
            if mul(r, mul(delta, delta, N), N):
                continue
            f = t ^ mul(alpha, delta, N)
            u = t ^ delta
            G = graph(n, g, f)
            C = graph(n, r, u)
            assert G ^ C == (1 << (n * (n - 1) // 2)) - 1
            assert G & 1
            fixed.add(G)
            raw += 1
        if r == 0:
            expected = 1
        else:
            k = (r & -r).bit_length() - 1
            d = N - k
            expected = 1 << (d // 2)
        assert len(fixed) == expected, (n, r, len(fixed), expected)
        per_r[r] = len(fixed)
        result.update(fixed)
    by_d = defaultdict(int)
    for r in per_r:
        if r:
            k = (r & -r).bit_length() - 1
            by_d[N - k] += 1
    assert all(by_d[d] == 1 << (d - 1) for d in range(1, N))
    expected_total = 1 + sum(1 << (d - 1 + d // 2) for d in range(1, N))
    assert len(result) == expected_total
    return result, raw


def sum_formula(n):
    assert n >= 1
    return 1 + sum(1 << (3 * k // 2) for k in range(n - 1))


def a113836(m):
    return sum(1 << (i + i // 2) for i in range(m))


def main():
    print('Source: Mwafise, OEIS A401496 revision 9; definition A401493 revision 10.')
    print('Fidelity: reproduce total labeled Riordan count (4^(n-1)+2)/3,')
    print('OEIS complement counts, and Pascal adjacency from binomial parity.')
    print('n presentations distinct_Riordan complement_compatible unit_side canonical_presentations')
    for n in range(1, MAX_ENUM + 1):
        N = n - 1
        graphs = exhaustive(n)
        full = (1 << (n * (n - 1) // 2)) - 1
        compatible = {G for G in graphs if full ^ G in graphs}
        total = (4 ** N + 2) // 3
        assert 3 * total == 4 ** N + 2
        assert len(graphs) == total, (n, len(graphs), total)
        assert len(compatible) == OEIS[n - 1] == sum_formula(n)
        if n == 1:
            units, raw = {0}, 1
        else:
            units, raw = classified_units(n)
            assert units == {G for G in compatible if G & 1}
            assert compatible == units | {full ^ G for G in units}
            assert not units & {full ^ G for G in units}
        presentations = 1 if n == 1 else (1 << N) * (1 << (N - 1))
        print(n, presentations, len(graphs), len(compatible), len(units), raw)
        H = (1 << N) - 1
        pascal = graph(n, H, mul(2, H, N)) if n > 1 else 0
        expected_pascal = sum((comb(i - 1, j) % 2) << (i * (i - 1) // 2 + j)
                              for i in range(1, n) for j in range(i))
        assert pascal == expected_pascal

    assert tuple(a113836(m) for m in range(1, 14)) == A113836
    assert tuple(sum_formula(n) for n in range(1, 15)) == OEIS
    # Obtain coefficients independently by division of the source polynomials:
    # (1-x)(1-8x^2) = 1-x-8x^2+8x^3, numerator x+x^2-6x^3.
    numerator = {1: 1, 2: 1, 3: -6}
    gf_coeff = [0] * 201
    for n in range(1, len(gf_coeff)):
        gf_coeff[n] = numerator.get(n, 0) + sum(
            c * gf_coeff[n - j] for j, c in [(1, 1), (2, 8), (3, -8)] if j <= n)
        an = sum_formula(n)
        assert gf_coeff[n] == an
        assert sum_formula(n + 1) - an == 1 << (3 * (n - 1) // 2)
        assert an == a113836(n - 1) + 1
        if n >= 3:
            assert an == 8 * sum_formula(n - 2) - 4
        if n % 2:
            k = (n - 1) // 2
            assert an == Fraction(3 * 8 ** k + 4, 7)
        else:
            k = n // 2
            assert an == Fraction(10 * 8 ** (k - 1) + 4, 7)
    print('All 14 published A401496 values and 13 A113836 values agree.')
    print('Recurrence, difference, sum, A113836, intended parity forms, and g.f.: n=1..200 OK.')
    # Literal source typo is exposed rather than silently corrected.
    literal_odd = Fraction(38 ** 1 + 4, 7)
    literal_even = Fraction(108 ** 1 + 4, 7)
    assert literal_odd == 6 != OEIS[2]
    assert literal_even == 16 != OEIS[3]
    print('SOURCE_TYPO: literal 38^1 gives a(3)=6, but direct count is 4;')
    print('SOURCE_TYPO: literal 108^1 gives a(4)=16, but direct count is 12.')
    print('The intended products 3*8^k and 10*8^(k-1) pass; no literal-typo refutation is proposed.')
    print('Finite diagnostics only; no all-order proof or Lean probe was performed.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
