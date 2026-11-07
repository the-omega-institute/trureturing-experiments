#!/usr/bin/env python3
"""Exact checks of Mikhin--Xiourouppa's odd Gaussian-splitting coefficient.

Only integer and Fraction arithmetic is used. No floating-point approximation,
computer algebra, or Lean proof is used. The reported finite-spacing numerical
enclosures are rational interval computations with alternating Taylor bounds.
"""

from fractions import Fraction as Q
from math import comb, factorial


def solve(matrix, vector):
    """Exact Gaussian elimination over the rationals or rational intervals."""
    a = [list(row) + [b] for row, b in zip(matrix, vector)]
    n = len(a)
    for k in range(n):
        pivot = a[k][k]
        a[k] = [v / pivot for v in a[k]]
        for i in range(n):
            if i != k:
                multiplier = a[i][k]
                a[i] = [x - multiplier * y for x, y in zip(a[i], a[k])]
    return [row[-1] for row in a]


def gamma_half(n):
    """Gamma(n+1/2)/sqrt(pi)."""
    return Q(factorial(2 * n), 4**n * factorial(n))


def integral_i(n, u):
    """I_n(a), u = (1+a)^(-1/2), after t = u^(-2)-1."""
    return 2 * sum(
        (Q((-1)**k * comb(n, k), 2 * k + 3) * (1 - u**(2 * k + 3))
         for k in range(n + 1)), Q(0))


def closed_coefficient(m, sigma, u):
    """The proposed coefficient multiplied by sqrt(pi)."""
    return (-Q(m, 12) * gamma_half(m + 2) * integral_i(m - 1, u)
            * integral_i(m, u) / (sigma * factorial(m - 1)))


def source_coefficients(m, sigma, u):
    """Source (16)--(26): all Gram and b coefficients scaled by sqrt(pi).

    beta is (17), alpha is (21), G and right sides are (26),
    and the returned L0,L2 are (23),(24). The component index is 1..m.
    This is independent of the Laguerre/subspace-velocity calculation.
    """
    def beta(n):
        return Q((-1)**n, factorial(n) * 4**n) / sigma * (u**(2 * n + 1) - 1)

    def alpha(i, j):
        n = i + j
        return Q((-1)**n * comb(2 * n, 2 * i),
                 factorial(n) * 2**(2 * n - 1)) / sigma

    g = [[alpha(i, j) for j in range(1, m + 1)] for i in range(1, m + 1)]
    s0 = solve(g, [beta(i) for i in range(1, m + 1)])
    node_sum = sum(j * j for j in range(1, m + 1))
    exceptional = beta(m + 1) - sum(
        (alpha(m + 1, j) * s0[j - 1] for j in range(1, m + 1)), Q(0))
    r2 = [node_sum * ((exceptional if i == m else 0) - alpha(i, m + 1) * s0[-1])
          for i in range(1, m + 1)]
    s2 = solve(g, r2)
    f2 = Q(1, 2) - u / sigma + 1 / (2 * sigma)
    l0 = f2 - sum((beta(i) * s0[i - 1] for i in range(1, m + 1)), Q(0))
    l2 = -sum((beta(i) * s2[i - 1] for i in range(1, m + 1)), Q(0))
    l2 -= beta(m + 1) * s0[-1] * node_sum
    return l0, l2


def projection_limit(m, sigma, u):
    """Independent limiting polynomial Gram projection in weight x^-1/2 exp(-x)."""
    g = [[gamma_half(i + j) for j in range(1, m + 1)] for i in range(1, m + 1)]
    b = [gamma_half(i) * (u**(2 * i + 1) - 1) for i in range(1, m + 1)]
    c = solve(g, b)
    return (sigma - 2 * u + 1 - sum((x * y for x, y in zip(b, c)), Q(0))) / (2 * sigma)


def add_term(d, k, v):
    d[k] = d.get(k, 0) + v
    if d[k] == 0:
        del d[k]


def bilinear(left, right):
    """Expand source inner products; 'T' denotes the standard normal target."""
    result = {}
    for i, a in left.items():
        for j, b in right.items():
            if i == 'T' and j == 'T':
                key = ('target-square', 0)
            elif i == 'T' or j == 'T':
                key = ('d', abs(j if i == 'T' else i))
            else:
                key = ('c', abs(i - j))
            add_term(result, key, a * b)
    return result


def check_source_equation_10():
    residual = {'T': 1, 0: -1}
    assert bilinear(residual, residual) == {
        ('target-square', 0): 1, ('d', 0): -2, ('c', 0): 1}
    for m in range(1, 9):
        fm = {m: 1, 0: -2, -m: 1}
        assert bilinear(residual, fm) == {
            ('c', 0): 2, ('c', m): -2, ('d', 0): -2, ('d', m): 2}
        for k in range(1, 9):
            expected = {}
            for key, value in [(('c', abs(k - m)), 2), (('c', k + m), 2),
                               (('c', m), -4), (('c', k), -4), (('c', 0), 4)]:
                add_term(expected, key, value)
            assert bilinear(fm, {k: 1, 0: -2, -k: 1}) == expected
    print('FIDELITY equation (10): F-square, 8 b entries, 64 A entries agree exactly')


class Interval:
    """Closed rational interval; all operations round outwards exactly."""
    def __init__(self, lo, hi=None):
        self.lo, self.hi = Q(lo), Q(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def cast(x):
        return x if isinstance(x, Interval) else Interval(x)

    def __add__(self, other):
        b = self.cast(other)
        return Interval(self.lo + b.lo, self.hi + b.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) + -self

    def __mul__(self, other):
        b = self.cast(other)
        products = [x * y for x in [self.lo, self.hi] for y in [b.lo, b.hi]]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.cast(other)
        assert b.lo * b.hi > 0, 'interval denominator contains zero'
        return self * Interval(1 / b.hi, 1 / b.lo)

    def __rtruediv__(self, other):
        return self.cast(other) / self


def exp_neg(z):
    assert Q(0) <= z <= Q(1)
    upper = sum(((-z)**k / factorial(k) for k in range(41)), Q(0))
    return Interval(upper - z**41 / factorial(41), upper)


def source_minimum_interval(m, sigma, u, s):
    """sqrt(pi)*L² at positive s=(h/sigma)², directly from source (9)--(11)."""
    def c(i, j):
        return exp_neg(Q((i - j)**2, 4) * s) / (2 * sigma)

    def d(i):
        return (u / (2 * sigma)) * exp_neg(i**2 * sigma**2 * s / (2 * (1 + sigma**2)))

    a = [[2 * (c(k, j) + c(-k, j) - 2 * c(0, j) - 2 * c(0, k) + 2 * c(0, 0))
          for k in range(1, m + 1)] for j in range(1, m + 1)]
    b = [2 * (c(0, 0) - c(0, j) - d(0) + d(j)) for j in range(1, m + 1)]
    weights = solve(a, b)
    return Q(1, 2) - 2 * d(0) + c(0, 0) - sum(
        (w * v for w, v in zip(weights, b)), Interval(0))


def display_interval(value, digits=20):
    scale = 10**digits
    lo = value.lo.numerator * scale // value.lo.denominator
    hi = -((-value.hi.numerator * scale) // value.hi.denominator)
    return f'[{lo}/{scale}, {hi}/{scale}]'


def main():
    check_source_equation_10()
    parameters = [(Q(1, 7), Q(1, 5)), (Q(7, 17), Q(7, 13))]
    for sigma, u in parameters:
        assert u**2 == 2 * sigma**2 / (1 + sigma**2)
        for m in range(1, 9):
            l0, l2 = source_coefficients(m, sigma, u)
            assert l0 == projection_limit(m, sigma, u)
            assert l2 == closed_coefficient(m, sigma, u)
            assert l2 < 0
            print(f'EXACT M={m} sigma={sigma}: sqrt(pi)*L2={l2}; source=closed; negative')
    for m in [1, 2]:
        for sigma, u in parameters:
            s = Q(1, 10**6)
            l0 = projection_limit(m, sigma, u)
            coarse = (source_minimum_interval(m, sigma, u, s) - l0) / s
            fine = (source_minimum_interval(m, sigma, u, s / 4) - l0) / (s / 4)
            richardson = (4 * fine - coarse) / 3
            predicted = closed_coefficient(m, sigma, u)
            tolerance = abs(predicted) / 10**8
            assert predicted - tolerance < richardson.lo <= richardson.hi < predicted + tolerance
            assert coarse.hi < 0 and fine.hi < 0
            print(f'RATIONAL_NUMERICAL M={m} sigma={sigma} s=1/1000000,1/4000000')
            print(f'  difference quotient at first s: {display_interval(coarse)}')
            print(f'  difference quotient at second s: {display_interval(fine)}')
            print(f'  Richardson: {display_interval(richardson)}; relative discrepancy < 1/100000000')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
