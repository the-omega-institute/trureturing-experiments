#!/usr/bin/env python3
"""Exact checks for Aivazidis, arXiv:2606.03320v1, Question 10.1.

Polynomial tuples store coefficients in increasing degree. No floating-point
arithmetic or external packages are used. Finite checks do not prove the
unbounded assertion. Run with Python 3; no input files are required.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import prod


def add(a, b):
    return tuple((a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b))))


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return tuple(c)


def shift(a, d):
    assert d >= 0
    return (0,) * d + a


def evaluate(a, x):
    z = 0
    for c in reversed(a):
        z = z * x + c
    return z


@lru_cache(None)
def gaussian(n, k):
    """Gaussian polynomial by the q-Pascal recurrence."""
    if k < 0 or k > n:
        return (0,)
    if k == 0 or k == n:
        return (1,)
    return add(gaussian(n - 1, k), shift(gaussian(n - 1, k - 1), n - k))


def partitions(n, ceiling=None):
    if n == 0:
        yield ()
        return
    for h in range(min(n, n if ceiling is None else ceiling), 0, -1):
        for tail in partitions(n - h, h):
            yield (h,) + tail


def conjugate(lam):
    return tuple(sum(x >= j for x in lam) for j in range(1, max(lam, default=0) + 1))


def widths(lam):
    return {h: lam[h - 1] - (lam[h] if h < len(lam) else 0)
            for h in range(1, len(lam) + 1)}


def parameters(lam):
    w = widths(lam)
    eta = sum(h * h // 4 for h in conjugate(lam))
    theta = prod(d + 1 for h, d in w.items() if h % 2)
    return w, eta, theta


def contained_columns(heights, i=0, ceiling=None):
    if i == len(heights):
        yield ()
        return
    top = heights[i] if ceiling is None else min(heights[i], ceiling)
    for x in range(top + 1):
        for tail in contained_columns(heights, i + 1, x):
            yield (x,) + tail


def count_polynomial(lam):
    """Proposition 2.1, summed over all subgroup column types."""
    heights = conjugate(lam)
    total = (0,)
    for xs in contained_columns(heights):
        term = (1,)
        for i, (h, x) in enumerate(zip(heights, xs)):
            y = xs[i + 1] if i + 1 < len(xs) else 0
            term = mul(term, shift(gaussian(h - y, x - y), y * (h - x)))
        total = add(total, term)
    return total


def paired_polynomial(lam):
    """Proposed normalized polynomial lower bound F_lam(q)."""
    w = widths(lam)
    ans = (w.get(1, 0) + 1,)
    for a in range(1, (len(lam) + 1) // 2 + 1):
        e, d = w.get(2 * a + 1, 0), w.get(2 * a, 0)
        odd = tuple(e for _ in range(a + 1))
        even = mul((1,) * a, (0,) + (1,) * d)
        ans = mul(ans, add(add((1,), odd), even))
    return ans


def beta(p, h):
    if h < 2:
        return Fraction(1)
    q = Fraction(1, p)
    c = sum(q ** j for j in range(h // 2 + 1))
    return 2 * c / (1 + c) if h % 2 == 0 else (1 + c) / 2


def paired_selected_types(lam):
    """Enumerate the proposed subfamily, independently of all-type recursion."""
    w = widths(lam)
    pairs = [(a, w.get(2 * a + 1, 0), w.get(2 * a, 0))
             for a in range((len(lam) + 1) // 2, 0, -1)]
    for ks in product(*(range(e + d + 1) for a, e, d in pairs)):
        prefix = tuple(x for (a, e, d), k in zip(pairs, ks)
                       for x in ((a + 1,) * k + (a,) * (e + d - k)))
        for t in range(w.get(1, 0) + 1):
            yield prefix + (1,) * t + (0,) * (w.get(1, 0) - t)


def direct_subgroup_count(lam, p):
    """Count actual subgroups by adjoining elements and taking cyclic cosets.

    This uses only addition in the finite product of residue class groups,
    and does not use partitions, Gaussian coefficients or Delsarte's formula.
    """
    moduli = tuple(p ** a for a in lam)
    elems = tuple(product(*(range(m) for m in moduli)))
    index = {x: i for i, x in enumerate(elems)}
    table = tuple(tuple(index[tuple((a + b) % m for a, b, m in zip(x, y, moduli))]
                        for y in elems) for x in elems)
    zero = index[(0,) * len(lam)]
    start = frozenset([zero])
    seen, queue = {start}, [start]
    for H in queue:
        for x in range(len(elems)):
            if x in H:
                continue
            K, y = set(H), x
            while y not in H:
                K.update(table[h][y] for h in H)
                y = table[y][x]
            K = frozenset(K)
            if K not in seen:
                seen.add(K)
                queue.append(K)
    return len(seen)


def main():
    source_lam = (5, 5, 3, 1)
    expected = (15, 13, 24, 33, 38, 41, 40, 37, 30, 14, 3)
    assert count_polynomial(source_lam) == expected
    assert parameters(source_lam)[1:] == (10, 3)
    print('FIDELITY: source Section 2 example lambda=(5,5,3,1)')
    print('s(A_lambda) coefficients in increasing degree:', expected)

    small = [lam for n in range(13) for lam in partitions(n)]
    assert len(small) == 272
    coefficient_cases = selected_cases = rational_cases = 0
    for lam in small + [source_lam]:
        w, eta, theta = parameters(lam)
        poly = count_polynomial(lam)
        assert len(poly) == eta + 1 and poly[-1] == theta
        normalized = tuple(reversed(poly))
        lower = paired_polynomial(lam)
        gap = add(normalized, tuple(-x for x in lower))
        assert all(x >= 0 for x in gap), (lam, gap)
        assert lower[0] == theta
        if len(lam) >= 2:
            assert len(lower) >= 2 and lower[1] > 0
        coefficient_cases += 1
        heights = conjugate(lam)
        selected = list(paired_selected_types(lam))
        assert len(selected) == len(set(selected))
        for xs in selected:
            assert len(xs) == len(heights)
            assert all(0 <= x <= h for x, h in zip(xs, heights))
            assert all(xs[i] >= xs[i + 1] for i in range(len(xs) - 1))
            selected_cases += 1
        for p in (2, 3, 5, 7):
            q = Fraction(1, p)
            exact = evaluate(poly, p)
            paired = p ** eta * evaluate(lower, q)
            factors = prod(beta(p, h) for h, d in w.items() if h >= 2 and d > 0)
            height_bound = theta * p ** eta * factors
            assert all(beta(p, h) > 1 for h in range(2, len(lam) + 1))
            assert exact >= paired >= height_bound >= theta * p ** eta
            if len(lam) < 2:
                assert exact == paired == height_bound
            else:
                assert height_bound > theta * p ** eta
            rational_cases += 1
    print('COEFFICIENTWISE: %d partitions (all sizes <=12, plus source example)' % coefficient_cases)
    print('PAIRED FAMILY: %d distinct valid selected types checked' % selected_cases)
    print('RATIONAL BOUNDS: %d cases at primes 2,3,5,7' % rational_cases)

    direct_cases = 0
    for p, max_size in ((2, 5), (3, 3), (5, 2)):
        for n in range(max_size + 1):
            for lam in partitions(n):
                direct = direct_subgroup_count(lam, p)
                assert direct == evaluate(count_polynomial(lam), p), (p, lam, direct)
                direct_cases += 1
    print('INDEPENDENT GROUP ENUMERATION: %d cases; |A|<=32 (p=2), <=27 (p=3), <=25 (p=5)' % direct_cases)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
