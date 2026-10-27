#!/usr/bin/env python3
"""Exact finite checks for Gaunt's Conjecture 2.4; no quadrature or floats.

The analytic half-integer formula is classical (DLMF 10.39.2, 10.49.12):
K_(m+1/2)(x) = sqrt(pi/(2*x))*exp(-x)*P_m(1/x).
We check coefficient identities, not the validity of that analytic formula.
"""
from fractions import Fraction as Q
from math import factorial


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(p, q):
    return trim([(p[i] if i < len(p) else Q(0))
                 + (q[i] if i < len(q) else Q(0))
                 for i in range(max(len(p), len(q)))])


def scale(c, p):
    return trim([c * a for a in p])


def shift(p, k):
    return [Q(0)] * k + p


def derivative(p):
    return trim([i * p[i] for i in range(1, len(p))] or [Q(0)])


def half_poly(m):
    return [Q(factorial(m + k), factorial(k) * factorial(m - k) * 2**k)
            for k in range(m + 1)]


def negative_endpoint_derivative(nu, m):
    # Divide -d[x^-nu K_(m+1/2)(x)]/dx by
    # sqrt(pi/2)*exp(-x)*x^(-nu-1/2), with z=1/x.
    p = half_poly(m)
    return add(add(p, scale(nu + Q(1, 2), shift(p, 1))),
               shift(derivative(p), 2))


def main():
    assert not any(isinstance(x, float) for x in half_poly(20))

    # Fidelity: source Appendix B recurrence (B.3), independently comparing
    # the factorial formula with the three-term recurrence.
    for m in range(1, 21):
        assert half_poly(m + 1) == add(
            half_poly(m - 1), scale(Q(2*m + 1), shift(half_poly(m), 1)))
    print('FIDELITY: Appendix B K recurrence, half-integer orders: 20 exact identities')

    # Fidelity: source (B.5), both signs of the order, including nu=-1/2.
    # The polynomial index of K_(nu+1) uses K_-r=K_r.
    for sign in [-1, 1]:
        for m in range(21):
            nu = sign * Q(2*m + 1, 2)
            j = abs(nu + 1) - Q(1, 2)
            assert j.denominator == 1 and j >= 0
            assert negative_endpoint_derivative(nu, m) == half_poly(int(j))
    print('FIDELITY/EQUALITY: Appendix B weighted derivative and n=1 identity: 42 exact identities')

    # Expected settlement: all pairs of half-integer a,b with A<=12 and
    # 0<=B<=A-2. The difference is a polynomial with zero constant/linear
    # terms and positive remaining coefficients: it is positive for EVERY
    # x>0 in these 66 admissible parameter cases, not only at sample points.
    cases = 0
    for A in range(2, 13):
        for B in range(A - 1):
            a, b = Q(2*A + 1, 2), Q(2*B + 1, 2)
            nu = (a*a - b*b - 1) / 2
            n = a - nu
            d = 1 - n
            assert d > 0 and nu - d >= 0 and (nu-d)**2 >= 2*d
            assert a*a - 2*nu - 1 == b*b
            gap = add(half_poly(A), scale(-1, negative_endpoint_derivative(nu, B)))
            assert gap[:2] == [0, 0]
            assert len(gap) >= 3 and all(c > 0 for c in gap[2:])
            cases += 1
    assert cases == 66
    print('SETTLEMENT: strict derivative gap for every x>0 in 66 half-integer parameter pairs')
    print('EXAMPLE: (n,nu,a,b)=(0,5/2,5/2,1/2); normalized gap=3/x^2')

    # Includes discriminant zero, negative n, and positive n. Squaring the
    # source threshold is legitimate because nu-d>=0 is checked separately.
    geom = 0
    for q in [Q(1,4), Q(1,2), Q(1), Q(2), Q(3)]:
        for h in [Q(0), Q(1,4), Q(1), Q(4)]:
            d, n = 2*q*q, 1-2*q*q
            nu = d + 2*q + h
            a = nu + n
            D = a*a - 2*nu - 1
            assert nu > 0 and a > 1 and D >= 0
            assert nu*nu - D == d*(2*nu+n+1) > 0
            assert (a-1)**2 - D == 2*d > 0
            if h == 0:
                assert D == 0
            geom += 1
    print(f'GEOMETRY: {geom} rational boundary/interior cases, including b=0')

    # Secant/tangent gaps of monomials used in the positive power series.
    # These finite checks do not certify an infinite sum or an integral.
    secants = tangents = 0
    for b in [Q(1,4), Q(1,2), Q(1), Q(2), Q(4)]:
        for delta in [Q(1,4), Q(1), Q(3)]:
            nu = b + delta
            aa, bb, cc = (b-1)**2, (b+1)**2, b*b+2*nu+1
            assert 0 <= aa < bb < cc
            for j in range(2, 21):
                gap = cc**j-bb**j-(cc-bb)*(bb**j-aa**j)/(bb-aa)
                assert gap > 0
                secants += 1
    for nu in [Q(1,4), Q(1), Q(2), Q(5)]:
        C = 2*nu+1
        for j in range(2, 21):
            assert C**j - 1 - j*(C-1) > 0
            tangents += 1
    print(f'CONVEXITY: {secants} secant and {tangents} b=0 tangent monomial checks')

    # Source Remark 2.5: the first two normalized tail terms agree with K_b;
    # the third tail coefficient exceeds the endpoint coefficient.
    asymptotic = 0
    for A in range(2, 13):
        for B in range(A-1):
            a, b = Q(2*A+1,2), Q(2*B+1,2)
            nu = (a*a-b*b-1)/2
            n = a-nu
            k1 = lambda r: (4*r*r-1)/8
            k2 = lambda r: (4*r*r-1)*(4*r*r-9)/128
            assert k1(a)-nu-Q(1,2) == k1(b)
            tail2 = k2(a)-(nu+Q(3,2))*k1(a)+(nu+Q(1,2))*(nu+Q(3,2))
            assert tail2-k2(b) == (nu*nu-b*b)/2 == (1-n)*(2*nu+n+1)/2 > 0
            for c in [b+Q(1,4), b+1]:
                assert k1(c)-k1(b) == (c*c-b*b)/2 > 0
            asymptotic += 1
    print(f'SHARPNESS: {asymptotic} exact two-term matches, positive third gaps, and larger-order signs')
    print('LIMIT: finite exact algebra; no general real-parameter theorem or analytic remainder bound certified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
