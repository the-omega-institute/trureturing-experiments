#!/usr/bin/env python3
"""Exact finite checks for Han–Molnár Lemma 11's general C*-algebra question.

Uses integers, Fractions and SymPy exact matrices only. This is evidence for a
preregistration outline, not a proof of its infinite-dimensional conclusion.
"""
from fractions import Fraction as F
from itertools import product
from math import isqrt
import sympy as s


def equal(a, b):
    assert (a - b).applyfunc(s.simplify) == s.zeros(a.rows, a.cols)


def projection(v):
    return v * (v.conjugate().T * v).inv() * v.conjugate().T


def fidelity():
    # Lemma 11, printed pp.20–21: B #_(1-p) P = (Tr(B^-1 P))^-p P.
    count = 0
    for n in (2, 3):
        for entries in product((1, 2), repeat=n):
            root = s.diag(*entries)
            b = root.inv() ** 2
            for coords in product((-1, 0, 1), repeat=n):
                if not any(coords):
                    continue
                pmat = projection(s.Matrix(coords))
                r = root * pmat * root
                alpha = s.trace(root ** 2 * pmat)
                for p in (s.Rational(1, 2), s.Rational(1, 3), s.Rational(2, 3)):
                    # Matrix powers are independently evaluated by SymPy's exact spectral algebra.
                    lhs = root.inv() * (r ** (1 - p)) * root.inv()
                    equal(lhs, alpha ** (-p) * pmat)
                    count += 1
    print(f"FIDELITY rank-one formula: {count} exact cases")


def corner_checks():
    count = 0
    projections = [s.eye(2), s.zeros(2)]
    projections += [projection(s.Matrix(v)) for v in product((-1, 0, 1), repeat=2) if any(v)]
    projections += [projection(s.Matrix([1, s.I])), projection(s.Matrix([2, s.I]))]
    for a, b, c in product(range(1, 4), range(1, 4), (-1, 0, 1)):
        if a * b <= c * c:
            continue
        root = s.Matrix([[a, c], [c, b]])
        cmat = root ** 2
        for pmat in projections:
            # Invert the compression as a full matrix after adjoining identity on its complement.
            x = (pmat * cmat * pmat + s.eye(2) - pmat).inv() * pmat
            equal(x * cmat * x, x)
            equal(pmat * x, x)
            q = root * x * root
            equal(q * q, q)
            equal(q.conjugate().T, q)
            for p in (s.Rational(1, 2), s.Rational(1, 3), s.Rational(2, 3)):
                equal(root.inv() * (q ** (1 - p)) * root.inv(), x)
            # Positivity of the logarithmic Jensen integrand via a Schur complement.
            if s.trace(pmat) == 1:
                w = pmat.columnspace()[0]
                z = (s.eye(2) - pmat).columnspace()[0]
                wn, zn = (w.conjugate().T * w)[0], (z.conjugate().T * z)[0]
                for t in (s.Rational(0), s.Rational(1), s.Rational(3)):
                    ct = cmat + t * s.eye(2)
                    aa = (w.conjugate().T * ct * w)[0] / wn
                    ff = (z.conjugate().T * ct * z)[0] / zn
                    ee2 = s.simplify(abs((w.conjugate().T * ct * z)[0]) ** 2 / (wn * zn))
                    compressed_inverse = (w.conjugate().T * ct.inv() * w)[0] / wn
                    gap = s.simplify(compressed_inverse - 1 / aa)
                    assert gap == s.simplify(ee2 / (aa * (aa * ff - ee2)))
                    assert gap.is_nonnegative is True
                    assert (gap == 0) == (ee2 == 0)
            count += 1
    # Noncommuting projections with rank two in dimension three.
    v = s.Matrix([1, 1, 1])
    pmat = s.eye(3) - projection(v)
    for entries in product((1, 2), repeat=3):
        root = s.diag(*entries)
        cmat = root ** 2
        x = (pmat * cmat * pmat + s.eye(3) - pmat).inv() * pmat
        equal(x * cmat * x, x)
        q = root * x * root
        equal(q * q, q)
        equal(root.inv() * (q ** s.Rational(1, 2)) * root.inv(), x)
        count += 1
    print(f"SETTLEMENT corner-inverse and projection identities: {count} exact cases")
    print("STRICTNESS Schur-resolvent gap: exact rational nonnegativity and zero iff off-diagonal zero")


class Interval:
    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = F(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, other):
        other = box(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -box(other)

    def __rsub__(self, other):
        return box(other) + -self

    def __mul__(self, other):
        other = box(other)
        endpoints = [x * y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return Interval(min(endpoints), max(endpoints))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = box(other)
        assert other.lo > 0 or other.hi < 0
        return self * Interval(1 / other.hi, 1 / other.lo)


def box(x):
    return x if isinstance(x, Interval) else Interval(x)


def log_bounds(x):
    x = F(x)
    assert x > 0
    k = 0
    while x >= 2:
        x /= 2
        k += 1
    while x < 1:
        x *= 2
        k -= 1

    def near_one(y):
        z = (y - 1) / (y + 1)
        n = 40
        total = sum((2 * z ** (2 * j + 1) / (2 * j + 1) for j in range(n)), F(0))
        # Positive atanh series and an exact geometric majorant of its remaining tail.
        tail = 2 * z ** (2 * n + 1) / ((2 * n + 1) * (1 - z * z))
        return Interval(total, total + tail)
    return near_one(x) + k * near_one(F(2))


def sqrt_bounds(x):
    x = box(x)
    assert x.lo >= 0
    scale = 10 ** 12
    lo = isqrt(x.lo.numerator * scale ** 2 // x.lo.denominator)
    hi = isqrt(x.hi.numerator * scale ** 2 // x.hi.denominator) + 1
    assert F(lo, scale) ** 2 <= x.lo
    assert F(hi, scale) ** 2 >= x.hi
    return Interval(F(lo, scale), F(hi, scale))


def invertible_witness():
    # Original invertible-domain check: p=1/2, B=diag(1,1/4),
    # P=(1/2)[[1,1],[1,1]], X=(2/5)P+epsilon I.
    eps = F(1, 10000)
    trace_x = F(2, 5) + 2 * eps
    determinant_r = 4 * eps * (F(2, 5) + eps)
    d = sqrt_bounds(determinant_r)
    trace_mean = (trace_x + F(5, 4) * d) / sqrt_bounds(1 + 5 * eps + 2 * d)
    entropy = eps * log_bounds(eps) + (F(2, 5) + eps) * log_bounds(F(2, 5) + eps)
    entropy += (F(1, 5) + eps) * log_bounds(4)
    h = entropy - 2 * (trace_x - trace_mean)
    # Coarse rational enclosure keeps the output readable and proves a strict sign.
    assert F(-1, 10) < h.lo <= h.hi < F(-1, 20)
    print("INVERTIBLE witness: -1/10 < H(X) < -1/20 (exact rational interval certificate)")
    print("WITNESS p=1/2; B=diag(1,1/4); X=(1/5)[[1,1],[1,1]]+(1/10000)I")
    # Singular limit's logarithmic-compression gap: K=log(5/4)P > 0.
    assert log_bounds(F(5, 4)).lo > 0
    print("LIMIT H((2/5)P)=-(2/5)log(5/4)<0 (exact sign)")


def main():
    fidelity()
    corner_checks()
    invertible_witness()
    print("SCOPE finite matrix evidence; no proof of the general C*-algebra implication")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
