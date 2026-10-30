#!/usr/bin/env python3
"""Exact checks for Mathar's OEIS A174193 recurrence; no Lean proof is claimed."""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial
import sys
import sympy as sp

LIMIT = 200
LISTED = [
    1, 1, 2, 6, 23, 103, 488, 2320, 10867, 49971, 225861, 1005821,
    4423807, 19255967, 83093364, 355956164, 1515437339, 6417727259,
    27054820649, 113602953089, 475370846639, 1983142721039,
    8250982614809, 34246497698633, 141838788939623, 586322071478503,
    2419474311617803, 9968270142800955, 41010408788528703, 168499207575017087,
]


def coefficient(v, n):
    return v[n] if 0 <= n < len(v) else Fraction(0)


def multiply(a, b, size):
    return [sum((coefficient(a, k) * coefficient(b, n-k)
                 for k in range(n+1)), Fraction(0)) for n in range(size)]


def source_coefficients(limit):
    # Generalized binomial expansion of the branch (1-4x)^(1/2) with s(0)=1.
    s = [Fraction(1)]
    for k in range(1, limit+3):
        s.append(s[-1] * Fraction(-4*(3-2*k), 2*k))
    assert s == [Fraction(1)] + [Fraction(-2*comb(2*k-2, k-1), k)
                                 for k in range(1, limit+3)]
    assert multiply(s, s, limit+3) == [1, -4] + [0]*(limit+1)
    P = [1, -8, 13, 24, -48]
    B = [1, -6, 1, 34, -26, -4]
    bs = multiply(B, s, limit+3)
    numerator = [coefficient(P, n)-bs[n] for n in range(limit+3)]
    assert numerator[:3] == [0, 0, 2]
    # Cancel x^2, then divide by the unit 2(1-x)(1-4x)^2, constant term 2.
    denominator = [2, -18, 48, -32]
    a = []
    for n in range(limit+1):
        a.append((numerator[n+2] - sum(denominator[k]*coefficient(a, n-k)
                                      for k in range(1, 4))) / denominator[0])
    assert all(v.denominator == 1 for v in a)
    return [int(v) for v in a]


def hook_length_sum(n):
    # The explicit formula proved in Section 2 of arXiv:2007.15111v1.
    ans = Fraction(comb(2*n, n), n+1)
    for k in range(n//2, n-1):
        f = Fraction(factorial(n)*(2*k-n+2),
                     (k+2)*factorial(k)*(n-k)*factorial(n-k-2))
        assert f.denominator == 1
        ans += f*f
    assert ans.denominator == 1
    return int(ans)


def avoids321(p):
    return all(not(p[i] > p[j] > p[k]) for i, j, k in combinations(range(len(p)), 3))


def almost321_count(n):
    return sum(avoids321(p) or any(avoids321(p[:i]+p[i+1:]) for i in range(n))
               for p in permutations(range(n)))


def certificate_checks():
    x, t = sp.symbols('x t')
    p = 1-4*x
    P = 1-8*x+13*x**2+24*x**3-48*x**4
    B = 1-6*x+x**2+34*x**3-26*x**4-4*x**5
    D = 2*x**2*(1-x)*p**2
    q = [-(t+2)*(404489*t-1866368),
         3147815*t**2-8680082*t-14098754,
         2*(-2336439*t**2+652022*t+14514849),
         8*(-1044770*t**2+10585023*t-20370553),
         8*(642982*t-2004177)*(2*t-9)]
    U = (x-1)*(4*x-1)*(2571928*x**2+1125370*x-404489)
    V = 2*(1970080*x**4+17265612*x**3-8693734*x**2-1192226*x+528695)
    W = -4542008*x**4+15852688*x**3+12946274*x**2-19631021*x+3732736
    H = 4*(677079*x**3+634434*x**2-3811346*x+933184)
    operator = sp.Poly(sum(x**j*q[j].subs(t, t+j) for j in range(5)), t)
    assert all(sp.Poly(operator.nth(k)-v, x).is_zero
               for k, v in enumerate([W, V, U]))
    assert operator.degree() == 2

    def E(f):
        return sp.diff(f, x)*D-f*sp.diff(D, x)

    def K(f):
        return x*E(f)*D+x**2*(sp.diff(E(f), x)*D-2*E(f)*sp.diff(D, x))

    rational = U*K(P)+V*x*E(P)*D+W*P*D**2-H*D**3
    radical = U*p**2*K(B)+(V*p**2-4*x*p*U)*x*E(B)*D + (
        W*p**2-2*x*p*V-2*x*(1-2*x)*U)*B*D**2
    assert sp.Poly(rational, x).is_zero
    assert sp.Poly(radical, x).is_zero
    assert sp.degree(H, x) == 3
    print('Euler operator conversion: exact identity')
    print('Rational polynomial certificate: 0')
    print('Square-root polynomial certificate: 0')
    print('Inhomogeneous term H coefficients:', [int(H.expand().coeff(x, k)) for k in range(4)])


def main():
    print('Python:', sys.version.split()[0], 'SymPy:', sp.__version__)
    a = source_coefficients(LIMIT)
    assert a[:len(LISTED)] == LISTED
    print('OEIS source fidelity: all 30 displayed coefficients agree')
    assert all(a[n] == hook_length_sum(n) for n in range(26))
    print('Section 2 hook-length sum: agrees for n=0,...,25')
    direct = [almost321_count(n) for n in range(8)]
    assert direct == a[:8]
    print('Independent permutation counts n=0,...,7:', direct)
    print('First 12 coefficients:', a[:12])
    residuals = []
    for n in range(4, LIMIT+1):
        residual = (-(n+2)*(404489*n-1866368)*a[n]
                    +(3147815*n*n-8680082*n-14098754)*a[n-1]
                    +2*(-2336439*n*n+652022*n+14514849)*a[n-2]
                    +8*(-1044770*n*n+10585023*n-20370553)*a[n-3]
                    +8*(642982*n-2004177)*(2*n-9)*a[n-4])
        residuals.append(residual)
    assert residuals == [0]*(LIMIT-3)
    print('Conjectured recurrence: zero residual for n=4,...,200 (197 cases)')
    certificate_checks()
    print('Finite checks are corroboration; the symbolic certificates support the universal argument.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
