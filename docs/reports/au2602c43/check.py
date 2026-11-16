#!/usr/bin/env python3
"""Exact finite checks of Au, arXiv:2602.08721v4, Conjecture 4.3.

No finite rank check proves the all-degree statement. Only Python's standard
library is used; every mathematical operation is rational or integer arithmetic.
"""
from fractions import Fraction as Q
from math import comb, factorial, gcd, lcm

DEGREE = 8
LAST_N = 79
MONOMIALS = [(i, j, d-i-j) for d in range(DEGREE+1)
             for i in range(d+1) for j in range(d-i+1)]
INDEX = {m: i for i, m in enumerate(MONOMIALS)}
PREV = [[INDEX.get(tuple(v-(axis == k) for k, v in enumerate(m)))
         for axis in range(3)] for m in MONOMIALS]
BY_DEGREE = [[i for i, m in enumerate(MONOMIALS) if sum(m) == d]
             for d in range(DEGREE+1)]


def affine_mul(p, constant, slopes):
    return [constant*p[i] + sum((s*p[j] for s, j in zip(slopes, PREV[i])
                                if s and j is not None), Q(0))
            for i in range(len(p))]


def affine_div(p, constant, slopes):
    assert constant != 0
    out = []
    for i in range(len(p)):
        out.append((p[i] - sum((s*out[j] for s, j in zip(slopes, PREV[i])
                                if s and j is not None), Q(0))) / constant)
    return out


def product_step(p, n, transformed=False):
    # Pochhammer quotient at n+1 divided by the quotient at n.
    numerator = [(Q(1,2)+n, (-1,-3,0)), (Q(1,2)+n, (-1,3,0)),
                 (Q(1,2)+2*n, (1,3,3)), (Q(3,2)+2*n, (1,3,3))]
    denominator = [(Q(1)+2*n, (-2,0,0)), (Q(2)+2*n, (-2,0,0)),
                   (Q(1)+n, (0,0,1)), (Q(1)+n, (0,2,1))]
    if transformed:
        # x=a, y=b+c, z=b, hence (a,b,c)=(x,z,y-z).
        numerator = [(u, (v[0], v[2], v[1]-v[2])) for u,v in numerator]
        denominator = [(u, (v[0], v[2], v[1]-v[2])) for u,v in denominator]
    for constant, slopes in numerator:
        p = affine_mul(p, constant, slopes)
    for constant, slopes in denominator:
        p = affine_div(p, constant, slopes)
    return [v/9 for v in p]


def transformed_coefficients(p):
    """Independent binomial substitution c=y-z in a,b,c Taylor coefficients."""
    out = [Q(0)] * len(p)
    for (i,j,k), value in zip(MONOMIALS,p):
        for h in range(k+1):
            out[INDEX[i,h,j+k-h]] += value*comb(k,h)*(-1)**(k-h)
    return out


def primitive(row):
    g = gcd(*row)
    return [x//g for x in row] if g else row


class RationalRowSpace:
    """Fraction-free Gaussian row reduction with exact integer gcd reduction."""
    def __init__(self):
        self.rows = {}

    def add(self, row):
        denominator = lcm(*(v.denominator for v in row))
        v = primitive([int(x*denominator) for x in row])
        for pivot, basis in sorted(self.rows.items()):
            if v[pivot]:
                g = gcd(v[pivot], basis[pivot])
                a, b = basis[pivot]//g, v[pivot]//g
                v = primitive([a*x-b*y for x,y in zip(v,basis)])
        pivot = next((i for i,x in enumerate(v) if x), None)
        if pivot is not None:
            self.rows[pivot] = v

    def rank(self):
        return len(self.rows)


def pochhammer(q, n):
    out = Q(1)
    for j in range(n):
        out *= q+j
    return out


def main():
    spaces = [[RationalRowSpace() for _ in range(DEGREE+1)] for _ in range(2)]
    p = [Q(1)] + [Q(0)]*(len(MONOMIALS)-1)
    pt = p[:]
    for n in range(LAST_N+1):
        raw = affine_mul(p, Q(-8*n-1), (2,-6,-6))
        invariant = affine_mul(pt, Q(-8*n-1), (2,-6,0))
        # Source fidelity: coefficient zero equals minus the stated a_n p_0.
        a_n = Q(1,9)**n*pochhammer(Q(1,2),n)*pochhammer(Q(1,4),n) \
              *pochhammer(Q(3,4),n)/Q(factorial(n)**3)
        assert raw[0] == -(1+8*n)*a_n
        if n == 0:
            assert raw == affine_mul([Q(1)]+[Q(0)]*(len(p)-1), Q(-1), (2,-6,-6))
        assert invariant == transformed_coefficients(raw)
        assert all(v == 0 for (i,j,k),v in zip(MONOMIALS,invariant) if k % 2)
        for degree, cols in enumerate(BY_DEGREE):
            spaces[n//40][degree].add([raw[i] for i in cols])
        if n < LAST_N:
            p = product_step(p,n)
            pt = product_step(pt,n,transformed=True)
        if n in (0,39,79):
            print('completed exact sequence index', n, flush=True)
    expected = [(n//2+1)*((n+1)//2+1) for n in range(DEGREE+1)]
    # Fidelity includes every displayed coefficient of the conjectured series.
    assert expected[:6] == [1,2,4,6,9,12]
    for window, family in enumerate(spaces):
        ranks = [s.rank() for s in family]
        assert ranks == expected, (window,ranks,expected)
        print(f'exact Q-ranks, n={40*window}..{40*window+39}: {ranks}')
    counts = [sum(n-2*k+1 for k in range(n//2+1)) for n in range(DEGREE+1)]
    assert counts == expected
    print('source fidelity: S(0,0,0;n)=-(1+8n)a_n for n=0..79')
    print('source fidelity: displayed Hilbert coefficients 1,2,4,6,9,12 reproduced')
    print('reflection: substitution (a,b,c)=(x,z,y-z), odd z coefficients vanish')
    print('checked: parameter degrees 0..8; all 80 indices; two disjoint rank windows')
    print('scope: finite checks only; no all-degree proof is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
