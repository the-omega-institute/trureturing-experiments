#!/usr/bin/env python3
"""Exact finite checks for the Lewis--Souza proportional-window limit route.

Only integer and Fraction arithmetic is used. The general convergence assertion
is not proved by these finite checks; its proposed argument is in the draft.
"""

from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations


def window(alpha, n):
    return tuple(a for a in range(1, n + 1) if alpha * n < a)


def primes_to(k):
    return tuple(p for p in range(2, k + 1)
                 if all(p % d for d in range(2, p)))


def prime_coordinates(a, primes):
    exponents = []
    large = a
    for p in primes:
        e = 0
        while large % p == 0:
            large //= p
            e += 1
        exponents.append(e)
    return (large, -large, *exponents)


@lru_cache(maxsize=None)
def dimension(elements):
    """Minimum realizer size, by all linear extensions and exact set cover.

    For each ordered incomparable pair (a,b), a realizer must contain an
    extension putting a after b. Comparabilities hold in every extension.
    Extension masks therefore give exactly a finite minimum set-cover problem.
    The zero-fold product convention gives dimension zero for size <= 1.
    """
    size = len(elements)
    if size <= 1:
        return 0
    pairs = tuple((i, j) for i in range(size) for j in range(size)
                  if elements[j] % elements[i] and elements[i] % elements[j])
    if not pairs:
        return 1
    if len(pairs) == size * (size - 1):
        # The increasing and decreasing natural ranks form an explicit
        # realizer; any incomparable pair excludes a one-chain embedding.
        coordinates = tuple((i, -i) for i in range(size))
        assert all((elements[j] % elements[i] == 0) ==
                   all(x <= y for x, y in zip(coordinates[i], coordinates[j]))
                   for i in range(size) for j in range(size))
        return 2
    predecessors = tuple(sum(1 << j for j in range(size)
                             if j != i and elements[i] % elements[j] == 0)
                         for i in range(size))
    masks = set()
    rank = [0] * size

    def extensions(used, level):
        if level == size:
            masks.add(sum(1 << k for k, (i, j) in enumerate(pairs)
                          if rank[i] > rank[j]))
            return
        for i in range(size):
            if not (used >> i & 1) and predecessors[i] & used == predecessors[i]:
                rank[i] = level
                extensions(used | (1 << i), level + 1)

    extensions(0, 0)
    full = (1 << len(pairs)) - 1
    assert masks and full == reduce(int.__or__, masks)
    if any(full ^ m in masks for m in masks):
        return 2
    by_bit = tuple(tuple(m for m in masks if m >> bit & 1)
                   for bit in range(len(pairs)))

    @lru_cache(maxsize=None)
    def can_cover(remaining, budget):
        if not remaining:
            return True
        if not budget:
            return False
        if any(remaining & m == remaining for m in masks):
            return True
        if budget == 1:
            return False
        bit = min((k for k in range(len(pairs)) if remaining >> k & 1),
                  key=lambda k: len(by_bit[k]))
        options = sorted({m & remaining for m in by_bit[bit]},
                         key=int.bit_count, reverse=True)
        return any(can_cover(remaining & ~m, budget - 1) for m in options)

    for d in range(2, size + 1):
        if can_cover(full, d):
            return d
    raise AssertionError('Every finite poset has a finite realizer')


def main():
    # Source fidelity: Section 6 says alpha > 1/2 gives an antichain of
    # dimension 2. Its dimension statement concerns at least two vertices.
    fidelity_cases = 0
    for alpha in (Fraction(3, 5), Fraction(2, 3), Fraction(3, 4), Fraction(9, 10)):
        for n in range(1, 81):
            xs = window(alpha, n)
            assert all(b % a and a % b for a, b in combinations(xs, 2))
            if 2 <= len(xs) <= 8:
                assert dimension(xs) == 2
                fidelity_cases += 1
    # A second source fidelity result: the positive divisors of m have
    # dimension omega(m), for these nonsingleton divisor sets (Introduction).
    divisor_values = []
    for m in (2, 4, 6, 12, 30):
        ds = tuple(d for d in range(1, m + 1) if m % d == 0)
        omega = sum(m % p == 0 for p in primes_to(m))
        actual = dimension(ds)
        assert actual == omega
        divisor_values.append((m, actual))
    print('FIDELITY antichain_cases=', fidelity_cases, sep='')
    print('FIDELITY divisor_dimensions=', divisor_values, sep='')

    alphas = sorted({Fraction(p, q) for q in range(2, 13)
                     for p in range(1, q)})
    bounded_pairs = 0
    for alpha in alphas:
        # Exact ceiling of 1/alpha, with positive numerator.
        k = (alpha.denominator + alpha.numerator - 1) // alpha.numerator
        primes = primes_to(k)
        assert 1 / alpha <= k
        for n in range(1, 61):
            xs = window(alpha, n)
            coords = {a: prime_coordinates(a, primes) for a in xs}
            assert len(set(coords.values())) == len(xs)
            for a in xs:
                for b in xs:
                    assert (b % a == 0) == all(x <= y for x, y in zip(coords[a], coords[b]))
                    bounded_pairs += 1
    print('BOUND alphas=', len(alphas), ' n=1..60 ordered_pairs=', bounded_pairs, sep='')

    persistence_windows = 0
    persistence_pairs = 0
    max_threshold = 0
    for alpha in alphas:
        for n0 in range(1, 31):
            xs = window(alpha, n0)
            m0 = (alpha.numerator * n0) // alpha.denominator + 1
            delta = m0 - alpha * n0
            ratio = alpha * n0 / delta
            q0 = ratio.numerator // ratio.denominator + 1
            threshold = q0 * n0
            assert delta > 0 and q0 >= 1 and q0 * delta > alpha * n0
            max_threshold = max(max_threshold, threshold)
            for q in sorted({q0, q0 + 1, 2 * q0 + 1}):
                for s in range(n0):
                    n = q * n0 + s
                    assert n >= threshold and n // n0 == q
                    images = tuple(q * a for a in xs)
                    assert len(set(images)) == len(xs)
                    for a, image in zip(xs, images):
                        assert image - alpha * n >= q * delta - alpha * s > 0
                        assert alpha * n < image <= n
                        assert image > 0
                    for a in xs:
                        for b in xs:
                            assert (b % a == 0) == ((q * b) % (q * a) == 0)
                            persistence_pairs += 1
                    persistence_windows += 1
    print('PERSISTENCE n0=1..30 all residues s=0..n0-1 at q=Q,Q+1,2Q+1',
          ' windows=', persistence_windows, ' ordered_pairs=', persistence_pairs,
          ' max_N0=', max_threshold, sep='')

    # Exact dimensions for every small case in this declared rational grid;
    # at the guaranteed persistence threshold verify dimension monotonicity
    # whenever both windows have at most eight elements.
    small_alphas = sorted({Fraction(p, q) for q in range(2, 7)
                           for p in range(1, q)})
    dimension_cases = 0
    dimension_persistence_cases = 0
    for alpha in small_alphas:
        values = []
        for n in range(1, 11):
            xs = window(alpha, n)
            if len(xs) > 8:
                continue
            d = dimension(xs)
            values.append((n, d))
            k = (alpha.denominator + alpha.numerator - 1) // alpha.numerator
            assert d <= len(primes_to(k)) + 2
            dimension_cases += 1
        print('SMALL_DIM alpha=', alpha, ' values=', values, sep='')
        for n0 in range(1, 6):
            delta = (alpha * n0).__floor__() + 1 - alpha * n0
            q0 = (alpha * n0 / delta).__floor__() + 1
            for n in range(q0 * n0, q0 * n0 + n0):
                before, after = window(alpha, n0), window(alpha, n)
                if max(len(before), len(after)) <= 8:
                    assert dimension(before) <= dimension(after)
                    dimension_persistence_cases += 1
    print('SMALL_DIM cases=', dimension_cases,
          ' persistence_comparisons=', dimension_persistence_cases, sep='')
    print('SCOPE finite exact checks; the general all-real-alpha theorem is not certified here')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
