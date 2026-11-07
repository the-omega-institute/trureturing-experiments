#!/usr/bin/env python3
"""Exact finite checks for Kimberling's A389142 urn-likelihood conjecture.

Source data: OEIS A389142 revision 11 and A111097 revision 35.
These checks are finite evidence, not a proof for all sample lengths.
"""
from fractions import Fraction as Q
from itertools import product
from math import comb, prod
import sympy as sp


def midsequence(n):
    return (Q(n * n) + Q(n * (n + 1), 2)) // 3


def proposed(n):
    return n * (n - 1) // 2 + (2 * n) // 3


def probability(n, k):
    assert n >= 1 and k >= 1
    if k < n:
        return Q(0)
    return Q(n * prod(range(k - n + 1, k + 1)), k ** (n + 1))


def event_count(n, k):
    return sum(len(set(w[:n])) == n and w[n] in w[:n]
               for w in product(range(k), repeat=n + 1))


def h(n, x):
    return (x - 1) ** (n + 1) - x ** n * (x - n)


def trunc(m, r, t):
    return sum(((-1) ** j * comb(m, j) * t ** j
                for j in range(min(m, r) + 1)), Q(0))


def main():
    a389142 = [0, 0, 2, 5, 8, 13, 19, 25, 33, 42, 51, 62, 74, 86,
               100, 115, 130, 147, 165, 183, 203, 224, 245, 268, 292,
               316, 342, 369, 396, 425, 455, 485, 517, 550, 583, 618,
               654, 690, 728, 767, 806, 847, 889, 931, 975, 1020, 1065,
               1112, 1160, 1208, 1258, 1309, 1360, 1413]
    assert [midsequence(n) for n in range(len(a389142))] == a389142
    print('FIDELITY A389142 revision 11: all 54 displayed terms reproduced')
    for n in range(1, 6):
        for k in range(1, 7):
            assert Q(event_count(n, k), k ** (n + 1)) == probability(n, k)
    assert probability(3, 5) == Q(36, 125)
    vals = {k: probability(3, k) for k in range(1, 101)}
    maximum = max(vals.values())
    assert [k for k, p in vals.items() if p == maximum] == [5]
    print('FIDELITY A111097 example: n=3, unique maximizer k=5 in 1..100; P=36/125')
    print('COMBINATORIAL BRIDGE: exhaustive draw counts n=1..5, k=1..6 agree with falling-factorial probability')
    for n in range(2, 61):
        target = midsequence(n)
        vals = {k: probability(n, k) for k in range(1, 2 * target + 11)}
        maximum = max(vals.values())
        assert [k for k, p in vals.items() if p == maximum] == [target]
    print('FINITE MAXIMUM: n=2..60, every k in 1..(2*M+10), unique maximum M')
    for n in range(2, 301):
        target = proposed(n)
        assert target >= n and target == midsequence(n) == (3*n*n+n+1)//6
        assert target == (n*n+n-2)//2 - (n-1)//3
        assert h(n, target) > 0 and h(n, target + 1) < 0
    print('BOUNDARY SIGNS: H_n(M)>0>H_n(M+1) for n=2..300, exact integers')
    grid = [Q(0), Q(1, 7), Q(1, 2), Q(6, 7), Q(1)]
    for m in range(41):
        for t in grid:
            assert trunc(m, 3, t) <= (1-t)**m <= trunc(m, 4, t)
            for r in (3, 4):
                coefficient = comb(m, r) if r <= m else 0
                assert trunc(m+1, r, t) == (1-t)*trunc(m, r, t) + (-1)**r*coefficient*t**(r+1)
    print('TRUNCATION: exact rational bounds and Pascal recurrence, m=0..40, t=0,1/7,1/2,6/7,1')
    m = sp.symbols('m')
    ns = [3*m, 3*m+1, 3*m+2]
    ms = [(9*m*m+m)/2, (9*m*m+7*m)/2, (9*m*m+13*m+4)/2]
    ls = [m*(m+1)/2, m*(7*m+5)/2, (m+1)*(2*m+1)]
    us = [-(117*m**4-2*m**3+67*m*m-6*m+8)/8,
          -m*(9*m**3+10*m*m+11*m+2)/8,
          -(63*m**4+166*m**3+193*m*m+106*m+24)/8]
    for n, target, lower, upper in zip(ns, ms, ls, us):
        A, B, C = n*(n+1)/2, n*(n+1)*(n-1)/6, n*(n+1)*(n-1)*(n-2)/24
        assert sp.expand(-target**2 + A*target - B - lower) == 0
        x = target + 1
        assert sp.expand(-x**3 + A*x**2 - B*x + C - upper) == 0
    print('SYMBOLIC CERTIFICATES: all six residue-class polynomial identities exact')
    for n in range(2, 21):
        for k in range(n, n+31):
            ratio = probability(n, k+1)/probability(n, k)
            assert ratio == Q(k**(n+1), (k+1)**n*(k+1-n))
            x = k+1
            e = Q(h(n, x), x**n)
            assert e == n-sum((Q(x-1, x)**j for j in range(n+1)), Q(0))
    print('RATIO AND GEOMETRIC IDENTITY: n=2..20, k=n..n+30, exact rationals')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
