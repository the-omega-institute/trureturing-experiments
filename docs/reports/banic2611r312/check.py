#!/usr/bin/env python3
"""Exact independent checks for arXiv:2611.14735v1, Remarks 3.12–3.13.

Uses the actual uniform multinomial law. No numerical approximations are evidence.
Run with Python >=3.12 and SymPy >=1.13.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import factorial, prod
import sys
import sympy as sp

assert sys.version_info >= (3, 12)


def compositions(total, bins):
    if bins == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, bins - 1):
                yield (first,) + rest


def cyclic(k):
    return sum(abs(k[(i + 1) % len(k)] - x) for i, x in enumerate(k))


def linear(k):
    return sum(abs(b - a) for a, b in zip(k, k[1:]))


@lru_cache(None)
def occupancy_hist(n, N):
    out = Counter()
    for k in compositions(N, n):
        out[k] = factorial(N) // prod(factorial(x) for x in k)
    assert sum(out.values()) == n ** N
    return out


def expectation(hist, f):
    return Q(sum(w * f(k) for k, w in hist.items()), sum(hist.values()))


def variance(hist, f):
    return expectation(hist, lambda k: f(k) ** 2) - expectation(hist, f) ** 2


def v(n, N):
    return variance(occupancy_hist(n, N), cyclic) / N


def F(n):
    return sum(Q(a, n ** j) for j, a in enumerate([12, -96, 464, -1056, 960, -400], 1))


def main():
    # Independent labelled histories reproduce every multinomial mass in this grid.
    for n in range(3, 7):
        for N in range(1, 6):
            labelled = Counter()
            for history in product(range(n), repeat=N):
                labelled[tuple(history.count(i) for i in range(n))] += 1
            assert labelled == occupancy_hist(n, N), (n, N)
    print('LAW_BRIDGE: labelled-history counts equal N!/product(k_i!) for 3<=n<=6, 1<=N<=5')

    # Fidelity: source Section 2.2's exact raw edge moments and Section 2.3's
    # E(d^c)=n/(n-1) E(d), plus Lemma 3.1's cyclic variance decomposition.
    for n in range(3, 8):
        for N in range(1, 7):
            h = occupancy_hist(n, N)
            edge = lambda k: k[1] - k[0]
            adj = lambda k: k[2] - k[1]
            assert expectation(h, edge) == 0
            assert variance(h, edge) == Q(2 * N, n)
            assert expectation(h, lambda k: edge(k) * adj(k)) == Q(-N, n)
            if n >= 4:
                disj = lambda k: k[3] - k[2]
                assert expectation(h, lambda k: edge(k) * disj(k)) == 0
            assert expectation(h, cyclic) == Q(n, n - 1) * expectation(h, linear)
            absedge = lambda k: abs(edge(k))
            ea = expectation(h, absedge)
            adjacent_cov = expectation(h, lambda k: abs(edge(k)) * abs(adj(k))) - ea ** 2
            rhs = n * variance(h, absedge) + 2 * n * adjacent_cov
            if n >= 4:
                disjoint_cov = expectation(h, lambda k: abs(edge(k)) * abs(k[3] - k[2])) - ea ** 2
                rhs += n * (n - 3) * disjoint_cov
            assert variance(h, cyclic) == rhs
    print('FIDELITY: exact edge moments, mean transfer and Lemma 3.1 hold for 3<=n<=7, 1<=N<=6')

    for n in range(3, 25):
        h = occupancy_hist(n, 4)
        overlap = lambda k: sum(min(x, k[(i + 1) % n]) for i, x in enumerate(k))
        assert all(cyclic(k) == 8 - 2 * overlap(k) for k in h)
        assert expectation(h, overlap) == Q(12, n) - Q(24, n ** 2) + Q(20, n ** 3)
        assert expectation(h, lambda k: overlap(k) ** 2) == Q(12, n) + Q(48, n ** 2) - Q(112, n ** 3)
        assert v(n, 4) == F(n)
    print('FAMILY: overlap identity, both moments and v(n,4)=F(n) verified for 3<=n<=24')

    x, t = sp.symbols('x t', positive=True)
    f = sum(sp.Rational(a) / x ** j for j, a in enumerate([12, -96, 464, -1056, 960, -400], 1))
    em = 12 / x - 24 / x ** 2 + 20 / x ** 3
    em2 = 12 / x + 48 / x ** 2 - 112 / x ** 3
    assert sp.cancel(em2 - em ** 2 - f) == 0
    p = t ** 5 + 4 * t ** 4 + 20 * t ** 3 + 144 * t ** 2 + 336 * t + 120
    assert sp.cancel((-sp.diff(f, x) * x ** 7 / 12).subs(x, t + 4) - p) == 0
    assert all(c > 0 for c in sp.Poly(p, t).all_coeffs())
    print('EXACT_ALGEBRA: F=E(M^2)-E(M)^2; F\'(x)=-12*P(x-4)/x^7; P coefficients =', sp.Poly(p, t).all_coeffs())

    expected = {3: Q(644, 729), 4: Q(247, 256), 5: Q(108, 125)}
    for n, ans in expected.items():
        assert v(n, 4) == ans
        print(f'v({n},4)={ans}')
    assert v(4, 4) - v(3, 4) == Q(15199, 186624) > 0
    assert v(4, 4) - v(5, 4) == Q(3227, 32000) > 0
    print('BIN_REVERSAL:', v(4, 4) - v(3, 4), '; RIGHT_DROP:', v(4, 4) - v(5, 4))

    types = defaultdict(Counter)
    aggregate = Counter()
    for k, w in occupancy_hist(4, 5).items():
        ty = tuple(sorted((a for a in k if a), reverse=True))
        types[ty][cyclic(k)] += w
        aggregate[cyclic(k)] += w
    assert dict(aggregate) == {2: 240, 4: 240, 6: 440, 8: 40, 10: 64}
    for ty in sorted(types, reverse=True):
        print('N5_TYPE', ty, dict(sorted(types[ty].items())))
    print('N5_DISTRIBUTION', dict(sorted(aggregate.items())), 'TOTAL', sum(aggregate.values()))
    h = occupancy_hist(4, 5)
    assert expectation(h, cyclic) == Q(315, 64)
    assert expectation(h, lambda k: cyclic(k) ** 2) == Q(925, 32)
    assert variance(h, cyclic) == Q(19175, 4096)
    assert v(4, 5) == Q(3835, 4096)
    assert v(4, 5) - v(4, 4) == Q(-117, 4096)
    print('v(4,5)=', v(4, 5), '; SAMPLE_REVERSAL:', v(4, 5) - v(4, 4))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
