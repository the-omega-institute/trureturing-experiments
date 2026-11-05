#!/usr/bin/env python3
"""Exact rational interval checks for Bouafia's support-pruning algorithm.

All evidence uses integers, Fraction and symbolic algebra. Intervals are closed;
every elementary operation rounds outward to a fixed rational grid. No binary
floating point is used. This is a reproducible arithmetic certificate, not a
Lean proof. Point indices in printed output are one-based.
"""
from fractions import Fraction as F
from math import factorial, isqrt
import sympy as sp

PRECISION = 40
S = 10 ** PRECISION
TERMS = 40


def floor(q):
    return q.numerator // q.denominator


def ceil(q):
    return -floor(-q)


class I:
    def __init__(self, lo, hi=None):
        lo = F(lo)
        hi = lo if hi is None else F(hi)
        assert lo <= hi
        self.lo, self.hi = floor(lo * S), ceil(hi * S)

    @classmethod
    def raw(cls, lo, hi):
        a = cls(0)
        a.lo, a.hi = lo, hi
        assert lo <= hi
        return a

    @staticmethod
    def lift(a):
        return a if isinstance(a, I) else I(a)

    def __add__(a, b):
        b = I.lift(b)
        return I.raw(a.lo + b.lo, a.hi + b.hi)

    __radd__ = __add__

    def __neg__(a):
        return I.raw(-a.hi, -a.lo)

    def __sub__(a, b):
        return a + -I.lift(b)

    def __rsub__(a, b):
        return I.lift(b) + -a

    def __mul__(a, b):
        b = I.lift(b)
        products = [a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi]
        return I.raw(min(products) // S, -((-max(products)) // S))

    __rmul__ = __mul__

    def __truediv__(a, b):
        b = I.lift(b)
        assert b.hi < 0 or 0 < b.lo, 'division interval contains zero'
        reciprocal = I(F(S, b.hi), F(S, b.lo))
        return a * reciprocal

    def __pow__(a, k):
        assert isinstance(k, int) and k >= 0
        v = I(1)
        for _ in range(k):
            v = v * a
        return v

    def square(a):
        if a.lo <= 0 <= a.hi:
            return I.raw(0, ceil(F(max(a.lo * a.lo, a.hi * a.hi), S)))
        return a * a

    def sqrt(a):
        assert a.lo >= 0
        lo = isqrt(a.lo * S)
        hi = isqrt(a.hi * S)
        if hi * hi < a.hi * S:
            hi += 1
        assert lo * lo <= a.lo * S and a.hi * S <= hi * hi
        return I.raw(lo, hi)

    def contains(a, q):
        return F(a.lo, S) <= q <= F(a.hi, S)

    def subset(a, lo, hi):
        return F(lo) <= F(a.lo, S) and F(a.hi, S) <= F(hi)

    def display(a, digits=18):
        scale = 10 ** digits
        lo = floor(F(a.lo * scale, S))
        hi = ceil(F(a.hi * scale, S))
        def dec(k):
            sign = '-' if k < 0 else ''
            k = abs(k)
            return sign + str(k // scale) + '.' + str(k % scale).zfill(digits)
        return '[' + dec(lo) + ', ' + dec(hi) + ']'


def exp_negative(d):
    """Enclose exp(-d), d >= 0, using exp(d/4) then reciprocation.

    For 0 <= x <= 1 and n > 0, sum_(k<n) x^k/k! <= exp(x),
    exp(x) <= that sum + x^n*(n+1)/(n!*n). The upper bound is
    Mathlib Real.exp_bound'. The lower bound follows from nonnegative terms.
    """
    x = d / 4
    assert 0 <= x.lo and x.hi <= S
    l, u = F(x.lo, S), F(x.hi, S)
    lower = sum((l ** k / factorial(k) for k in range(TERMS)), F(0))
    upper = sum((u ** k / factorial(k) for k in range(TERMS)), F(0))
    upper += u ** TERMS * (TERMS + 1) / (factorial(TERMS) * TERMS)
    e = I(lower, upper)
    return (I(1) / e) ** 4


def kernel(points):
    n = len(points)
    z = [[I(1) for _ in points] for _ in points]
    for i in range(n):
        for j in range(i):
            r = sum(((a - b).square() for a, b in zip(points[i], points[j])), I(0))
            z[i][j] = z[j][i] = exp_negative(r.sqrt())
    return z


def solve(z, active):
    """Interval Gaussian elimination; every pivot is proved strictly positive.

    Positive pivots of the symmetric principal matrix certify nonsingularity;
    every real matrix in the input box has its all-ones solve in the result.
    """
    a = [[z[i][j] for j in active] + [I(1)] for i in active]
    pivots = []
    n = len(active)
    for k in range(n):
        pivot = a[k][k]
        assert pivot.lo > 0, 'pivot not certified positive'
        pivots.append(pivot)
        for i in range(k + 1, n):
            factor = a[i][k] / pivot
            for j in range(k + 1, n + 1):
                a[i][j] = a[i][j] - factor * a[k][j]
            a[i][k] = I(0)
    v = [I(0) for _ in range(n)]
    for i in reversed(range(n)):
        v[i] = (a[i][n] - sum((a[i][j] * v[j] for j in range(i + 1, n)), I(0))) / a[i][i]
    w = [I(0) for _ in z]
    for i, value in zip(active, v):
        w[i] = value
    return w, pivots


def residual(z, w):
    return [sum((a * b for a, b in zip(row, w)), I(0)) - 1 for row in z]


COORDINATES = [('0.124', '0.511'), ('0.208', '0.314'), ('-0.923', '-0.097'),
               ('-0.762', '-0.517'), ('0.098', '-0.295'), ('0.286', '0.460'),
               ('-0.469', '-0.059'), ('0.918', '-0.319'), ('-0.349', '0.043')]


def check_case(label, coordinate_radius=F(0), scale_radius=F(0)):
    print(label)
    scale = I(1 - scale_radius, 1 + scale_radius)
    points = [[I(F(a) - coordinate_radius, F(a) + coordinate_radius) * scale
               for a in row] for row in COORDINATES]
    # Distinctness holds uniformly throughout the box.
    for i in range(9):
        for j in range(i):
            assert any(a.hi < b.lo or b.hi < a.lo for a, b in zip(points[i], points[j]))
    z = kernel(points)
    w0, pivots0 = solve(z, list(range(9)))
    assert all(w0[i].hi < 0 for i in [1, 6])
    assert all(w0[i].lo > 0 for i in [0, 2, 3, 4, 5, 7, 8])
    active = [i for i in range(9) if w0[i].lo > 0]
    assert active == [0, 2, 3, 4, 5, 7, 8]
    for i, w in enumerate(w0):
        print('first w_' + str(i + 1), w.display())
    print('retained', [i + 1 for i in active])
    w, pivots = solve(z, active)
    assert all(w[i].lo > 0 for i in active)
    assert w[1].lo == w[1].hi == w[6].lo == w[6].hi == 0
    for i in active:
        print('second w_' + str(i + 1), w[i].display())
    q = residual(z, w)
    assert q[1].lo > 0
    assert q[6].hi < floor(-F(14, 100000) * S)
    for i in active:
        assert q[i].contains(0)
    print('q_2', q[1].display())
    print('q_7', q[6].display())
    print('minimum pivot lower bound', min(p.lo for p in pivots0 + pivots), '/', S)
    mass = sum(w, I(0))
    assert mass.lo > 0
    # An explicit competitor, without relying on sufficiency of Proposition 3:
    # p = w/mass; v = (1-epsilon)*p + epsilon*delta_7 is a probability.
    # E(p)=1/mass and (Zp)_7 - E(p)=q_7/mass, using the exact solve.
    # E(v)-E(p) = 2e*q_7/mass + e^2*(1 - 1/mass - 2q_7/mass).
    epsilon = F(1, 100000)
    gap = 2 * epsilon * q[6] / mass + epsilon ** 2 * (1 - I(1) / mass - 2 * q[6] / mass)
    assert gap.hi < 0
    print('mass', mass.display())
    print('competitor energy difference', gap.display())
    # Literal printed guard stops at the first solve, returning negative masses.
    first_mass = sum(w0, I(0))
    assert first_mass.lo > 0
    assert (w0[1] / first_mass).hi < 0 and (w0[6] / first_mass).hi < 0
    print('printed-guard return has negative probabilities at 2 and 7')
    print('intended-guard return fails optimality at point 7')


def fidelity():
    a = sp.Symbol('a', real=True)
    z = sp.Matrix([[1, a], [a, 1]])
    inverse = sp.Matrix([[1, -a], [-a, 1]]) / (1 - a ** 2)
    assert sp.simplify(z * inverse - sp.eye(2)) == sp.zeros(2)
    w = sp.simplify(inverse * sp.ones(2, 1))
    assert all(sp.simplify(v - 1 / (1 + a)) == 0 for v in w)
    assert sp.simplify(sum(w) - (1 + (1 - a) / (1 + a))) == 0
    # The displayed two-point formula in Section 2.3 is reproduced identically.
    print('FIDELITY two-point inverse and magnitude 2/(1+a) = 1+(1-a)/(1+a): exact symbolic identity')
    print('Here a=exp(-L), L>0; (1-a)/(1+a)=tanh(L/2).')
    assert exp_negative(I(0)).contains(F(1))
    assert I(0).sqrt().contains(F(0))
    print('degenerate distances sqrt(0)=0 and exp(0)=1 enclosed')


def main():
    print('grid denominator', S, 'Taylor terms', TERMS)
    fidelity()
    check_case('EXACT NOMINAL NINE-POINT INSTANCE')
    check_case('UNIFORM PERTURBATION FAMILY: each input coordinate +/- 1/10^12, common scale 1 +/- 1/10^10',
               F(1, 10 ** 12), F(1, 10 ** 10))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
