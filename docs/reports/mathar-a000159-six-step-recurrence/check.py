#!/usr/bin/env python3
"""Exact independent checks of Mathar's first OEIS A000159 conjecture.

Sources: https://oeis.org/A000159/internal (revision 54), and
https://oeis.org/A058087/internal (revision 66). Python 3 and SymPy required.
The finite tests do not establish the universal recurrence. The rational
certificates check the algebra in the proposed formal-series argument.
"""

from fractions import Fraction
from math import comb, factorial
import sympy as sp


def rook_factor(n, k):
    assert n >= 1 and 0 <= k <= n
    q, r = divmod(2 * n * comb(2 * n - k, k), 2 * n - k)
    assert r == 0
    return q


def hit_coeff(n, degree):
    if n == 0:
        return int(degree == 0)
    return sum(
        (-1) ** (k - degree) * rook_factor(n, k)
        * factorial(n - k) * comb(k, degree)
        for k in range(degree, n + 1)
    )


def coefficients(n):
    return [
        2 * (-252307 * n + 1041077),
        504614 * n**2 - 3362985 * n + 5118150,
        1280831 * n**2 - 7397886 * n + 6461565,
        746598 * n**2 - 2913543 * n - 1336090,
        -405481 * n**2 + 6175011 * n - 15469320,
        -375862 * n**2 + 4098537 * n - 8846430,
        2 * (-187931 * n + 560630),
    ]


def assert_zero(expr):
    assert sp.cancel(expr) == 0


def main():
    # Fidelity: all six explicitly displayed A058087 polynomials, not a
    # sequence generated from the conjectured recurrence itself.
    rows = [[1], [-1, 2], [0, 0, 2], [1, 0, 3, 2],
            [2, 8, 4, 8, 2], [13, 30, 40, 20, 15, 2]]
    for n, row in enumerate(rows):
        assert [hit_coeff(n, r) for r in range(n + 1)] == row
    print('FIDELITY: A058087 displayed U_0,...,U_5 agree coefficientwise')

    a = [hit_coeff(n, 3) for n in range(301)]
    assert a[3:9] == [2, 8, 20, 152, 994, 7888]
    print('FIDELITY: A000159 a(3),...,a(8) =', a[3:9])
    for n in range(9, 301):
        assert sum(c * a[n - j] for j, c in enumerate(coefficients(n))) == 0
    print('RECURRENCE: all 292 source-domain indices 9 <= n <= 300 pass')
    for n in range(6, 9):
        assert sum(c * a[n - j] for j, c in enumerate(coefficients(n))) == 0
    print('ZERO EXTENSION: indices 6,7,8 also pass')

    # Separate factorial pullback expansion: B = sum k! x^k/(1+x)^(2k).
    bcoeff = [1] + [sum((-1) ** (n - k) * factorial(k)
                       * comb(n + k - 1, n - k)
                       for k in range(1, n + 1)) for n in range(1, 101)]
    x, z, d, n = sp.symbols('x z d n')
    rnum = -2 * (2*x**5 - 3*x**4 - 5*x**3 - 4*x**2 - 3*x - 2)
    snum = 2 * (2*x**5 - x**4 - 8*x**3 - 7*x**2 - 5*x - 2)
    denom = 3 * (1 + x)**2

    def conv_coeff(p, seq, i):
        p = sp.Poly(p, x)
        return sum(int(p.nth(j)) * seq[i-j]
                   for j in range(min(i, p.degree()) + 1))

    for i in range(101):
        assert conv_coeff(denom, a, i) == (
            conv_coeff(rnum, bcoeff, i) + int(sp.Poly(snum, x).nth(i)))
    print('PULLBACK: 3(1+x)^2 A = rnum B + snum through degree 100')

    for nn in range(1, 41):
        for k in range(nn + 1):
            j = nn - k
            term2 = comb(2*k + j - 1, j - 1) if j else 0
            assert rook_factor(nn, j) == comb(2*k + j, j) + term2
    print('REINDEXING: rook/binomial decomposition for 1 <= n <= 40 passes')

    # theta F = (1/z-1)F-1/z. Check the three displayed theta iterates.
    expected = [
        (1/z - 1, -1/z),
        (1 - 3/z + 1/z**2, 2/z - 1/z**2),
        (-1 + 7/z - 6/z**2 + 1/z**3, -3/z + 5/z**2 - 1/z**3),
    ]
    h, q = sp.Integer(1), sp.Integer(0)
    for eh, eq in expected:
        h, q = z*sp.diff(h, z) + h*(1/z - 1), -h/z + z*sp.diff(q, z)
        assert_zero(h - eh)
        assert_zero(q - eq)
    print('FACTORIAL SERIES: all three theta identities pass exact algebra')

    R, S = rnum / denom, snum / denom
    Z = x / (1+x)**2
    theta_combo = (4-4*x)*expected[2][0] + (18-6*x)*expected[1][0] + (20-2*x)*expected[0][0] + 6
    theta_const = (4-4*x)*expected[2][1] + (18-6*x)*expected[1][1] + (20-2*x)*expected[0][1]
    assert_zero(x**3/(3*(1+x)**4)*theta_combo.subs(z, Z) - R)
    assert_zero(x**3/(3*(1+x)**4)*theta_const.subs(z, Z) - S)
    bb = (1-x**3) / (x**2*(1+x))
    cc = (x**2-1) / x**2
    assert_zero(sp.diff(Z, x)*(1-Z)/Z**2 - bb)
    assert_zero(-sp.diff(Z, x)/Z**2 - cc)
    print('REDUCTION: diagonal-to-RB+S and chain-rule coefficients pass')

    P2 = -x**3*(x+1)*(375862*x**3+29619*x**2-776217*x-504614)
    P1 = -x*(375862*x**6+35945*x**5-2525682*x**4-2312643*x**3+993731*x**2+1849143*x+504614)
    P0 = -1133912*x**6+2249705*x**5+2743028*x**4-3357337*x**3-3210883*x**2+2259779*x+2082154
    rhs = 1136624*x**3-10348*x**4+1516100*x**5
    u, v = sp.diff(R, x) + R*bb, R*cc + sp.diff(S, x)
    assert_zero(P2*(sp.diff(u, x)+u*bb) + P1*u + P0*R)
    assert_zero(P2*(u*cc+sp.diff(v, x)) + P1*v + P0*S - rhs)
    print('CERTIFICATE: both rational differential identities pass')

    op = sp.Poly(sp.expand(sum(x**j * c.subs(n, d+j)
                              for j, c in enumerate(coefficients(n)))), d)
    assert_zero(x**2*op.nth(2) - P2)
    assert_zero(x*(op.nth(2) + op.nth(1)) - P1)
    assert_zero(op.nth(0) - P0)
    assert sp.degree(rhs, x) == 5
    print('OPERATOR: all three coefficients match the conjectured recurrence')
    # Fractions also check the coefficient extraction before theta rewriting.
    for k in range(41):
        assert Fraction(4*k**3+18*k**2+20*k+6, 3) == comb(2*k+3, 3) + comb(2*k+2, 2)
        assert Fraction(-4*k**3-6*k**2-2*k, 3) == -comb(2*k+3, 3) + comb(2*k+2, 2)
    print('COEFFICIENT EXTRACTION: theta weights for 0 <= k <= 40 pass')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
