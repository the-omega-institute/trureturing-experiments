#!/usr/bin/env python3
"""Exact finite checks for Svinin arXiv:2611.02309v1, Section 6.

No output of this program is an all-index proof. The direct coupled recurrence
is implemented independently of the ratio-map certificate and closed forms.
Optional gmpy2 accelerates exact integers; the default fallback is Python int.
"""
import argparse
from collections import deque
from fractions import Fraction as Q
from math import gcd as python_gcd
import sys
import sympy as sp

try:
    import gmpy2
    Integer = gmpy2.mpz
    gcd = gmpy2.gcd
    BACKEND = 'gmpy2 ' + gmpy2.version()
except ImportError:
    Integer, gcd, BACKEND = int, python_gcd, 'Python int'

A = (0, 0, 1, 0, 2, 1, 2, 3, 3, 5)
C = (0, 0, -1, 4, -1, 1, 1, 0, 4, 0)
W = (0, -1, 5, -5, 2, 0, -1, 4, -4, 2)
Z = (-1, 1, 0, -1, 2)


def valuations(n):
    k, r = divmod(n, 10)
    vx = 5*k*k + r*k + A[r]
    return vx, vx + 2*k + C[r]


def v2int(z):
    assert z != 0
    z = abs(z)
    return int((z & -z).bit_length() - 1)


def v2(q):
    q = Q(q)
    assert q
    return v2int(q.numerator) - v2int(q.denominator)


def modunit(q, modulus=8):
    q = Q(q)
    assert q.denominator % 2
    return q.numerator * pow(q.denominator, -1, modulus) % modulus


def F(a, b):
    return a*a*b*b + 8*a*a*b + 4*a*b*b - 39*a*b + 16*a + 8*b + 8


def next_ratio(x, y):
    return (4*y + 2)/(x*y*(y + 4))


def chart_check():
    charts = {
        'E': ((1, 3), (0, -1, 5, -5, 2, 0, -1),
              (0, -1, 2, -5, 4, 0, -1), (7, 6, 10, 0, 7, 3, 4), (3, 1)),
        'O': ((3, 1), (0, -1, 4, -4, 2, 0, -1),
              (0, -1, 2, -4, 3, 0, -1), (7, 6, 10, 2, 8, 5, 5), (1, 3)),
    }
    rows = 0
    for phase, (residue, vv, hh, pp, endpoint) in charts.items():
        pairs = [(a, b) for a in range(1, 128, 2) for b in range(1, 128, 2)
                 if (a % 8, b % 8) == residue and F(a, b) % 128 == 0]
        expected = ([(8*t+1, 8*t+3) for t in range(16)] if phase == 'E'
                    else [((9*(8*t+1)+2) % 128, 8*t+1) for t in range(16)])
        assert set(pairs) == set(expected) and len(pairs) == 16
        for a, b in sorted(pairs):
            orbit = [Q(a), Q(b, 2)]
            prec = [7, 6]
            for j in range(5):
                x, y = orbit[-2:]
                px, py = prec[-2:]
                vx, vy, h = v2(x), v2(y), v2(y+4)
                assert px > vx and py > max(vy, h)
                q = (4*y+2)/(y*(y+4))
                pq = py + min(v2(4*y*y+4*y+8), v2(4*y+2)+py) - 2*vy - 2*h
                pn = min(pq-vx, v2(q)+px-2*vx)
                orbit.append(next_ratio(x, y))
                prec.append(pn)
            assert tuple(map(v2, orbit)) == vv
            assert tuple(v2(r+4) for r in orbit) == hh
            assert tuple(prec) == pp
            assert all(p > max(v2(r), v2(r+4)) for p, r in zip(prec, orbit))
            assert (modunit(orbit[5]), modunit(2*orbit[6])) == endpoint
            rows += 1
        print('CHART', phase, 'classes', len(pairs), 'valuation', vv,
              'shifted', hh, 'precision', pp, 'endpoint', endpoint, flush=True)
    assert rows == 32
    # The derivative is odd for every odd pair, not just the selected roots.
    assert all((2*a*b*b+16*a*b+4*b*b-39*b+16) % 2 == 1
               for a in range(1, 128, 2) for b in range(1, 128, 2))
    print('CHART_CERTIFICATE_OK classes=32 modulus=128', flush=True)


def algebra_check():
    v, w, u, yp, x, xp, k = sp.symbols('v w u yp x xp k')
    inv = lambda v, w: 2*v*v*w*w+8*v*v*w+8*v*w*w-39*v*w+8*v+8*w+4
    t = (4*w+2)/(v*w*(w+4))
    assert sp.cancel(inv(w, t) - (4*w+2)*inv(v, w)/(v*v*w*(w+4))) == 0
    q = lambda w: (4*w+2)/(w*(w+4))
    assert sp.cancel(q(yp)-q(w) -
                     (w-yp)*(4*w*yp+2*(w+yp)+8)/(w*yp*(w+4)*(yp+4))) == 0
    assert inv(1, sp.Rational(3, 2)) == 0
    for r in range(10):
        def xx(j):
            shift, s = divmod(r+j, 10)
            kk = k+shift
            return 5*kk*kk+s*kk+A[s]
        def dd(j):
            shift, s = divmod(r+j, 10)
            return 2*(k+shift)+C[s]
        assert sp.expand(xx(3)+xx(0)-xx(1)-xx(2)-Z[r % 5]) == 0
        assert sp.expand(dd(1)-dd(0)-W[r]) == 0
    print('EXACT_ALGEBRA_OK invariant, error identity, all 10 symbolic residue recurrences', flush=True)


def rational_cross_check(limit=100):
    xs, ys = [Q(1), Q(1), Q(2)], [Q(1), Q(1), Q(3)]
    for n in range(limit-2):
        ys.append((4*xs[n+1]*ys[n+2]+2*xs[n+2]*ys[n+1])/xs[n])
        xs.append((4*ys[n+1]*xs[n+2]+ys[n+2]*xs[n+1])/ys[n])
    for n, (x, y) in enumerate(zip(xs, ys)):
        assert x > 0 and y > 0 and x.denominator == y.denominator == 1
        assert (v2(x), v2(y)) == valuations(n)
        assert python_gcd(x.numerator, y.numerator) == 2**min(valuations(n))
    # A separate scalar Somos-5 implementation, not the coupled update.
    for initial, direct in [([1, 1, 2, 11, 164], xs), ([1, 1, 3, 16, 194], ys)]:
        s = list(initial)
        for n in range(limit-4):
            quotient, rem = divmod(87*s[n+2]*s[n+3]+14*s[n+1]*s[n+4], s[n])
            assert rem == 0
            s.append(quotient)
        assert all(Q(a) == b for a, b in zip(s, direct))
    assert (4*Q(1)*Q(1)+Q(1)*Q(1))/Q(2) == Q(5, 2)  # y_(-1); backward integrality fails.
    print('RATIONAL_AND_SOMOS5_CROSS_CHECK_OK n=0..' + str(limit), flush=True)


def integer_orbit(limit):
    expected_x = [1,1,2,11,164,4210,274828,51195896,19837002616]
    expected_y = [1,1,3,16,194,6892,559512,80551288,34641433216]
    expected_g = [1,1,1,1,2,2,4,8,8]
    xs, ys = deque(map(Integer,[1,1,2]), maxlen=3), deque(map(Integer,[1,1,3]), maxlen=3)
    prev, initial_gcds = -1, []
    for n in range(limit+1):
        if n < 3:
            x, y = xs[n], ys[n]
        else:
            # Direct equations only. Neither predicted valuation nor ratio map is used here.
            yn, ry = divmod(4*xs[1]*ys[2]+2*xs[2]*ys[1], xs[0])
            xn, rx = divmod(4*ys[1]*xs[2]+ys[2]*xs[1], ys[0])
            assert rx == ry == 0, ('nonintegral', n)
            xs.append(xn); ys.append(yn)
            x, y = xn, yn
        assert x > 0 and y > 0
        vx, vy = v2int(x), v2int(y)
        assert (vx, vy) == valuations(n), ('valuation', n, vx, vy)
        exponent = min(vx, vy)
        assert gcd(x,y) == Integer(1) << exponent, ('gcd', n)
        assert exponent >= prev, ('monotonicity', n)
        prev = exponent
        if n < 9:
            assert (x,y,gcd(x,y)) == (expected_x[n],expected_y[n],expected_g[n])
            initial_gcds.append(int(gcd(x,y)))
        if n % 100 == 0 or n == limit:
            print('ORBIT', n, 'v2x', vx, 'v2y', vy, 'gcd_exponent', exponent,
                  'x_bits', x.bit_length(), 'y_bits', y.bit_length(), flush=True)
    print('SOURCE_FIDELITY_OK exact x,y,g tables n=0..8; det M=14=2*7', flush=True)
    print('SOURCE_GCD_TABLE', initial_gcds, flush=True)
    print('DIRECT_INTEGER_CHECK_OK n=0..' + str(limit) +
          ' exact divisions, positivity, both valuation formulas, entire gcd, exponent monotonicity', flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--limit', type=int, default=2000)
    args = p.parse_args()
    assert args.limit >= 8
    print('EXACT_INTEGER_BACKEND', BACKEND, flush=True)
    chart_check()
    algebra_check()
    rational_cross_check()
    integer_orbit(args.limit)
    print('ALL_CHECKS_OK', flush=True)

if __name__ == '__main__':
    main()
