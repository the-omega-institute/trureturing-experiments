#!/usr/bin/env python3
"""Exact checks for Han--Zhan--Zhou--Zhu, arXiv:2610.15716v1, Eq. (63).

All evidence uses integers, Fraction intervals, or exact symbolic algebra.
The finite checks do not prove the assertion for arbitrary distributions.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb
import sympy as sp

ZERO = F(0)
ONE = F(1)


def series_bounds(x):
    """Bounds for log(x), 1 <= x <= 2, by the positive atanh series."""
    assert ONE <= x <= 2
    z = (x - 1) / (x + 1)
    total = ZERO
    power = z
    n = 24
    for j in range(n):
        total += 2 * power / (2 * j + 1)
        power *= z * z
    remainder = 2 * power / ((2 * n + 1) * (1 - z * z))
    return total, total + remainder


LOG2 = series_bounds(F(2))


@lru_cache(maxsize=None)
def log_bounds(x):
    assert isinstance(x, F) and x > 0
    shift = 0
    while x < 1:
        x *= 2
        shift -= 1
    while x >= 2:
        x /= 2
        shift += 1
    # Enclose the argument on a dyadic grid before the series, keeping the
    # certificate arithmetic small even for long truncated geometric laws.
    grid = 2 ** 72
    floor = (x.numerator * grid) // x.denominator
    lower = F(floor, grid)
    upper = min(F(floor + 1, grid), F(2))
    lo = series_bounds(lower)[0]
    hi = series_bounds(upper)[1]
    if shift >= 0:
        return lo + shift * LOG2[0], hi + shift * LOG2[1]
    return lo + shift * LOG2[1], hi + shift * LOG2[0]


def entropy_bounds(p):
    lo = hi = ZERO
    for x in p:
        assert x >= 0
        if x:
            a, b = log_bounds(x)
            lo -= x * b
            hi -= x * a
    return lo, hi


def outward_text(interval, digits=12):
    """Outward decimal rendering using integer arithmetic only."""
    def render_integer(n):
        sign = '-' if n < 0 else ''
        n = abs(n)
        return sign + str(n // 10 ** digits) + '.' + str(n % 10 ** digits).zfill(digits)
    a, b = interval
    scale = 10 ** digits
    low = (a.numerator * scale) // a.denominator
    high = -((-b.numerator * scale) // b.denominator)
    return '[' + render_integer(low) + ', ' + render_integer(high) + ']'


def normalized(weights):
    total = sum(weights)
    assert total > 0
    return [F(x, total) for x in weights]


def log_concave(p):
    return all(p[i] * p[i] >= p[i - 1] * p[i + 1] for i in range(1, len(p) - 1))


def difference(p):
    n = len(p)
    return [sum((p[i] * p[i - d] for i in range(n) if 0 <= i - d < n), ZERO)
            for d in range(1 - n, n)]


def exact_structure_check(p):
    assert all(x > 0 for x in p) and sum(p) == 1 and log_concave(p)
    n = len(p)
    q = difference(p)
    assert sum(q) == 1 and q == q[::-1]
    pp = sorted(p, reverse=True)
    qq = sorted(q, reverse=True) + [ZERO]
    rr = sorted([x / 2 for x in p for _ in range(2)], reverse=True)
    for k in range(1, n + 1):
        # For each antidiagonal s, the half-endpoint central selector has
        # exactly the mass of the k largest entries of i -> p_i p_(s-i).
        central_total = ZERO
        for s in range(2 * n - 1):
            weights = [(i, p[i] * p[s - i]) for i in range(n) if 0 <= s - i < n]
            central = sum((v * (ONE if abs(2 * i - s) < k else
                                  F(1, 2) if abs(2 * i - s) == k else ZERO)
                           for i, v in weights), ZERO)
            top = sum(sorted((v for _, v in weights), reverse=True)[:k], ZERO)
            assert central == top, (p, s, k)
            central_total += central
        q0 = q[n - 1]
        band = q0 + 2 * sum(q[n:n + k - 1], ZERO)
        if k < n:
            band += q[n - 1 + k]
        assert band == central_total
        assert sum(pp[:k], ZERO) <= band <= sum(qq[:2 * k], ZERO)
    for m in range(2 * n + 1):
        assert sum(rr[:m], ZERO) <= sum(qq[:m], ZERO)
    return q


def entropy_check(p):
    q = exact_structure_check(p)
    hp = entropy_bounds(p)
    hq = entropy_bounds(q)
    delta = hq[0] - hp[1], hq[1] - hp[0]
    assert delta[1] < LOG2[0], (p, delta)
    return delta


def geometric_fidelity():
    a, z = sp.symbols('a z', positive=True)
    c = (1 - a) / (1 + a)
    # Exact generating functions of the one-sided geometric and its difference.
    pgf = (1 - a) / (1 - a * z)
    q_abs_pgf = c * (1 + 2 * a * z / (1 - a * z))
    assert sp.simplify(pgf.subs(z, 1) - 1) == 0
    assert sp.simplify(q_abs_pgf.subs(z, 1) - 1) == 0
    mean = sp.diff(pgf, z).subs(z, 1)
    abs_mean = sp.diff(q_abs_pgf, z).subs(z, 1)
    # H(p)=-log(1-a)-E[X]log(a), H(q)=-log(c)-E[|D|]log(a).
    L, M, A = sp.symbols('L M A')
    hp = -L - mean * A
    hq = -L + M - abs_mean * A
    increment = M - a * A / (1 + a)
    assert sp.simplify(hq - hp - increment) == 0
    # Check the coefficient formula by independent geometric-series summation.
    for av in [F(1, 3), F(1, 2), F(9, 10), F(99, 100)]:
        for d in range(9):
            convolution = (1 - av) ** 2 * av ** d / (1 - av ** 2)
            closed = (1 - av) * av ** d / (1 + av)
            assert convolution == closed
    print('FIDELITY: geometric difference q_d=(1-a)/(1+a)*a^|d| and source entropy identity: exact algebra OK')


def geometric_increment(a):
    l1, l2 = log_bounds(1 + a), log_bounds(a)
    c = a / (1 + a)
    return l1[0] - c * l2[1], l1[1] - c * l2[0]


def main():
    geometric_fidelity()
    count = 0
    nonmonotone = 0
    for n in range(1, 7):
        for w in product(range(1, 5), repeat=n):
            p = normalized(w)
            if not log_concave(p):
                continue
            entropy_check(p)
            count += 1
            if not (all(w[i] <= w[i + 1] for i in range(n - 1)) or
                    all(w[i] >= w[i + 1] for i in range(n - 1))):
                nonmonotone += 1
    print('FINITE: lengths 1..6, every integer weight 1..4; log-concave cases=' + str(count) +
          ', nonmonotone=' + str(nonmonotone) + '; antidiagonals, all majorization prefixes, entropy bounds OK')
    named = {
        'point mass': [ONE],
        'uniform length 9': normalized([1] * 9),
        'binomial n=8': normalized([comb(8, k) for k in range(9)]),
        'asymmetric nonmonotone': normalized([1, 3, 6, 8, 7, 5, 2]),
    }
    for a, n in [(F(1, 2), 12), (F(9, 10), 20), (F(9, 10), 80)]:
        named['truncated geometric a=' + str(a) + ' n=' + str(n)] = normalized([a ** k for k in range(n)])
    for name, p in named.items():
        print('ENTROPY ' + name + ': increment in ' + outward_text(entropy_check(p)))
    previous = None
    for a in [F(1, 3), F(1, 2), F(9, 10), F(99, 100), F(999, 1000)]:
        delta = geometric_increment(a)
        assert delta[1] < LOG2[0]
        if previous:
            assert previous[1] < delta[0]
        previous = delta
        print('INFINITE geometric a=' + str(a) + ': increment in ' + outward_text(delta) +
              '; gap to log(2) in ' + outward_text((LOG2[0] - delta[1], LOG2[1] - delta[0])))
    print('LOG(2) in ' + outward_text(LOG2))
    print('SCOPE: finite enumeration and exact interval checks; no universal theorem or equality classification certified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
