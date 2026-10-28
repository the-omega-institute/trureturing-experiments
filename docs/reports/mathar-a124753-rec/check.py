#!/usr/bin/env python3
"""Exact checks for the seven-term recurrence conjecture in OEIS A124753.

Source: https://oeis.org/search?q=id:A124753&fmt=text, revision #71.
Requires Python 3 and SymPy. No floating point arithmetic is used.
"""
from fractions import Fraction
from functools import reduce
from math import comb, factorial
from operator import mul

import sympy as sp


def a(n):
    assert isinstance(n, int) and n >= 0
    m, s = divmod(n, 3)
    return Fraction((s + 1) * comb(4 * m + s, m), 3 * m + s + 1)


def c(n):
    return [
        8019*n*(n-1)*(n+1),
        17496*n*(n-1)*(n-3),
        2592*(3*n-5)*(n-1)*(3*n-16),
        216*(-224*n**3+48*n**2+3926*n-6331),
        576*(-288*n**3+2448*n**2-6558*n+5443),
        768*(-288*n**3+3600*n**2-14878*n+20375),
        -8192*(4*n-23)*(2*n-11)*(4*n-21),
    ]


def product(xs):
    return reduce(mul, xs, 1)


def falling(z, r):
    assert r >= 0
    return product(z-h for h in range(r))


def rising(z, r):
    assert r >= 0
    return product(z+h for h in range(1, r+1))


# Integer coefficient arrays, increasing by degree; independent of SymPy.
def add(p, q):
    out = [0] * max(len(p), len(q))
    for i, v in enumerate(p):
        out[i] += v
    for i, v in enumerate(q):
        out[i] += v
    return out


def times(p, q):
    out = [0] * (len(p)+len(q)-1)
    for i, v in enumerate(p):
        for j, w in enumerate(q):
            out[i+j] += v*w
    return out


def scale(k, p):
    return [k*x for x in p]


def iprod(ps):
    return reduce(times, ps, [1])


def coefficient_polynomials(s):
    n = [s, 3]
    nm = lambda k: [s-k, 3]
    cube, square = iprod([n]*3), times(n, n)
    cubic = lambda b: scale(b[0], add(add(scale(b[1], cube), scale(b[2], square)),
                                             add(scale(b[3], n), [b[4]])))
    return [
        scale(8019, iprod([n, nm(1), [s+1, 3]])),
        scale(17496, iprod([n, nm(1), nm(3)])),
        scale(2592, iprod([[3*s-5, 9], nm(1), [3*s-16, 9]])),
        cubic([216, -224, 48, 3926, -6331]),
        cubic([576, -288, 2448, -6558, 5443]),
        cubic([768, -288, 3600, -14878, 20375]),
        scale(-8192, iprod([[4*s-23, 12], [2*s-11, 6], [4*s-21, 12]])),
    ]


def main():
    published = [1,1,1,1,2,3,4,9,15,22,52,91,140,340,612,969,2394,4389,
                 7084,17710,32890,53820,135720,254475,420732,1068012,2017356,
                 3362260,8579560,16301164,27343888,70068713,133767543,225568798,
                 580034052,1111731933,1882933364,4855986044,9338434700]
    assert [a(n) for n in range(len(published))] == published
    print('FIDELITY: all 39 published OEIS data terms, n=0..38, agree')
    for n in range(1, 151):
        rhs = sum((a(3*k)*a(n-1-3*k) for k in range((n-1)//3+1)), Fraction(0))
        assert a(n) == rhs, n
    print('FIDELITY: Manyama convolution from the entry agrees for n=1..150')
    for n in range(6, 301):
        assert sum((v*a(n-j) for j, v in enumerate(c(n))), Fraction(0)) == 0, n
    print('RECURRENCE: exact rational residuals vanish for n=6..300 (295 cases)')

    expected_pairs = [
        [(0,0),(-1,2),(-1,1),(-1,0),(-2,2),(-2,1),(-2,0)],
        [(0,1),(0,0),(-1,2),(-1,1),(-1,0),(-2,2),(-2,1)],
        [(0,2),(0,1),(0,0),(-1,2),(-1,1),(-1,0),(-2,2)],
    ]
    bridges = 0
    M = sp.Symbol('m')
    for s in range(3):
        pairs = [divmod(s-j, 3) for j in range(7)]
        assert pairs == expected_pairs[s]
        ps, ips = [], []
        for j, (t, u) in enumerate(pairs):
            length = 8-j+t
            assert -2 <= t <= 0 and 0 <= u <= 2 and length >= 0
            ps.append((u+1)*falling(M, -t)*rising(4*M+s-8, length)*falling(3*M+s+1, j))
            ips.append(scale(u+1, iprod(
                [[-h, 1] for h in range(-t)] +
                [[s-8+h, 4] for h in range(1, length+1)] +
                [[s+1-h, 3] for h in range(j)])))
            assert len(ips[-1]) == 9
            assert [int(sp.expand(ps[-1]).coeff(M, k)) for k in range(9)] == ips[-1]
            for m in range(2, 51):
                n = 3*m+s
                numerator, denominator = 4*m+s-j+t, 3*m+s-j+1
                assert min(n-j, numerator, m+t, denominator, 4*m+s-8) >= 0
                C = Fraction(factorial(4*m+s-8), factorial(m)*factorial(3*m+s+1))
                P = (u+1)*falling(m,-t)*rising(4*m+s-8,length)*falling(3*m+s+1,j)
                factorial_form = Fraction((u+1)*factorial(numerator), factorial(m+t)*factorial(denominator))
                assert a(n-j) == factorial_form == C*P, (s,j,m)
                bridges += 1
        sympy_cs = c(3*M+s)
        integer_cs = coefficient_polynomials(s)
        for sc, ic in zip(sympy_cs, integer_cs):
            assert [int(sp.expand(sc).coeff(M, k)) for k in range(4)] == ic
        h = sp.expand(sum(sc*p for sc, p in zip(sympy_cs, ps)))
        hc = [int(h.coeff(M, k)) for k in range(12)]
        ih = reduce(add, [times(ic, ip) for ic, ip in zip(integer_cs, ips)], [0])
        assert len(ih) == 12 and ih == hc == [0]*12 and h == 0
        print(f'POLYNOMIAL: s={s}, degree <=11, coefficients m^0..m^11={hc}; SymPy and integer convolution agree')
    print(f'FACTORIAL BRIDGE: all 21 shifts, m=2..50 ({bridges} exact cases), agree')
    print('SCOPE: finite sequence checks and exact polynomial certificates; no Lean proof is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
