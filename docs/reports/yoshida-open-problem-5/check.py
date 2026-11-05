"""Exact fidelity and refutation checks for Yoshida Open Problem 5.

All numerical evidence uses Python integers/Fractions or SymPy exact algebra.
"""

from fractions import Fraction

import sympy as sp


def ext(a, j):
    return Fraction(0) if j < 0 else a[j]


def turan(a, r, k):
    return a[k] * a[k] - ext(a, k - r) * ext(a, k + r)


def source_example_30():
    # Yoshida's Example 30: a_k = 2^(-k^2), with a_j=0 for j<0.
    a = [Fraction(1, 2 ** (k * k)) for k in range(32)]
    for r in range(1, 8):
        t = [turan(a, r, k) for k in range(24)]
        assert all(x > 0 for x in t)
        for k in range(1, 20):
            assert t[k] * t[k] >= 4 * ext(t, k - 1) * ext(t, k + 1)
    # The input identity in the example is exact at every interior index.
    for k in range(1, 24):
        assert a[k] * a[k] == 4 * a[k - 1] * a[k + 1]
    print("source Example 30: exact 4-rapid decrease checked for r=1..7")


def family(q, L, k):
    s = Fraction(q * q - 1)
    A = Fraction((q + 1) * (q + 1), 2)
    return L * (s / L) ** k * A ** (-(k * (k - 1) // 2))


def family_refutation():
    q = 7
    L = Fraction(49)  # L>s makes the coefficients strictly decreasing initially.
    s = Fraction(q * q - 1)
    A = Fraction((q + 1) * (q + 1), 2)
    a = [family(q, L, k) for k in range(16)]
    assert all(x > 0 for x in a)
    assert a[1] < a[0]

    # Exact adjacent-ratio identity and entire-series ratio.
    for k in range(1, 12):
        assert a[k] * a[k] == A * a[k - 1] * a[k + 1]
    for k in range(10):
        assert a[k + 1] / a[k] == (s / L) * A ** (-k)
    assert sp.limit((sp.Rational(s.numerator, s.denominator) /
                     sp.Rational(L.numerator, L.denominator)) *
                    sp.Rational(A.denominator, A.numerator) ** sp.Symbol("k"),
                    sp.Symbol("k"), sp.oo) == 0

    r = 1
    c = 1 - A ** (-r * r)
    b = [turan(a, r, k) for k in range(12)]
    assert all(x > 0 for x in b)
    # The exact output minimum is the k=r ratio A^2*c; the boundary on its
    # left is A^2/c and the tail ratios are A^2.
    output_ratio_at_r = b[r] * b[r] / (b[r - 1] * b[r + 1])
    assert output_ratio_at_r == A * A * c
    assert b[0] == a[0] * a[0]
    assert b[1] == c * a[1] * a[1]

    b1_lower = s * s * (1 - 1 / A)
    assert b[1] >= b1_lower
    u = sp.symbols("u", nonnegative=True, integer=True)
    qsym = u + 7
    ssym = qsym**2 - 1
    Asym = (qsym + 1)**2 / 2
    certificate = sp.expand(sp.cancel(4 * (ssym**2 * (1 - 1 / Asym) / 2 - Asym**2)))
    assert certificate == u**4 + 24*u**3 + 196*u**2 + 592*u + 368
    assert certificate.subs(u, 0) > 0

    # B(x)=max(4,(1+sqrt(1+x))^2/2) is strictly larger than x/2.
    # The exact lower bound gives b1/2 > A^2, hence B(b1)>A^2.
    assert b1_lower / 2 > A * A
    assert b[1] / 2 > A * A
    assert output_ratio_at_r < b[1] / 2
    print("counterexample q=7, L=49, r=1: output ratio", output_ratio_at_r,
          "< b1/2", b[1] / 2)
    print("certificate:", certificate)


source_example_30()
family_refutation()
print("ALL_CHECKS_OK")
