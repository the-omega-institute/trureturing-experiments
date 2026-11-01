#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A183241 and A182963 conjectures.

Uses Python's standard library only. Finite checks support a proposed general
argument; they do not prove the all-index integrality assertions.
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, factorial, prod


@lru_cache(None)
def partitions(n, largest=None):
    """Each integer partition once, as a decreasing tuple of positive parts."""
    if n == 0:
        return ((),)
    if largest is None:
        largest = n
    return tuple(
        (first,) + rest
        for first in range(min(n, largest), 0, -1)
        for rest in partitions(n - first, first)
    )


def multinomial(parts):
    numerator = factorial(sum(parts))
    denominator = prod(factorial(a) for a in parts)
    assert denominator > 0 and numerator % denominator == 0
    return numerator // denominator


def power_sums(r, bound):
    return [sum(multinomial(p) ** r for p in partitions(n))
            for n in range(bound + 1)]


def weighted_partition_product(r, bound):
    """The independently expanded OEIS product for b_n / (n!)^r."""
    f = [Fraction(1)] + [Fraction(0)] * bound
    for d in range(1, bound + 1):
        weight = Fraction(1, factorial(d) ** r)
        for n in range(d, bound + 1):
            f[n] += weight * f[n - d]
    return [f[n] * factorial(n) ** r for n in range(bound + 1)]


def source_recurrence(b):
    a = [Fraction(1)]
    for n in range(1, len(b)):
        a.append(sum((b[k] * a[n - k] for k in range(1, n + 1)),
                     Fraction(0)) / n)
    return a


def multiply(f, g, bound):
    return [sum((f[k] * g[n - k] for k in range(n + 1)), Fraction(0))
            for n in range(bound + 1)]


def formal_exponential(b):
    """Expand sum_{j=0}^N L^j/j!, independently of the recurrence."""
    bound = len(b) - 1
    logarithm = [Fraction(0)] + [Fraction(b[n], n)
                                for n in range(1, bound + 1)]
    term = [Fraction(1)] + [Fraction(0)] * bound
    result = term[:]
    for j in range(1, bound + 1):
        term = [v / j for v in multiply(term, logarithm, bound)]
        result = [a + v for a, v in zip(result, term)]
    return result


def primitive_ledger(b):
    c = [0] * len(b)
    for n in range(1, len(b)):
        residual = b[n] - sum(d * c[d] for d in range(1, n) if n % d == 0)
        assert residual >= 0 and residual % n == 0, (n, residual)
        c[n] = residual // n
        assert b[n] == sum(d * c[d] for d in range(1, n + 1) if n % d == 0)
    return c


def integral_euler_product(c):
    bound = len(c) - 1
    f = [1] + [0] * bound
    for d in range(1, bound + 1):
        factor = [0] * (bound + 1)
        factor[0] = 1
        if c[d]:
            for k in range(1, bound // d + 1):
                factor[d * k] = comb(c[d] + k - 1, k)
        f = multiply(f, factor, bound)
    return f


@lru_cache(None)
def content_words(content):
    if not any(content):
        return ((),)
    words = []
    for letter, count in enumerate(content):
        if count:
            smaller = list(content)
            smaller[letter] -= 1
            words.extend((letter,) + w for w in content_words(tuple(smaller)))
    return tuple(words)


def rotate(words, d):
    return tuple(w[d:] + w[:d] for w in words)


def least_period(words):
    n = len(words[0])
    return next(d for d in range(1, n + 1) if rotate(words, d) == words)


def word_spaces(r, bound):
    spaces = {}
    for n in range(1, bound + 1):
        elements = set()
        for content in partitions(n):
            words = content_words(content)
            assert len(words) == multinomial(content)
            elements.update((len(content), ws) for ws in product(words, repeat=r))
        assert len(elements) == power_sums(r, n)[n]
        spaces[n] = elements
    return spaces


def check_roots_and_orbits(r, bound, ledger):
    spaces = word_spaces(r, bound)
    primitives = {n: {x for x in xs if least_period(x[1]) == n}
                  for n, xs in spaces.items()}
    for n, xs in spaces.items():
        strata = {d: set() for d in range(1, n + 1) if n % d == 0}
        for l, words in xs:
            d = least_period(words)
            assert n % d == 0
            root = tuple(w[:d] for w in words)
            q = n // d
            assert tuple(w * q for w in root) == words
            content = tuple(Counter(root[0])[i] for i in range(l))
            assert all(k > 0 for k in content)
            assert content == tuple(sorted(content, reverse=True))
            assert all(tuple(Counter(w)[i] for i in range(l)) == content for w in root)
            assert (l, root) in primitives[d]
            strata[d].add((l, words))
        for d, stratum in strata.items():
            image = {(l, tuple(w * (n // d) for w in ws))
                     for l, ws in primitives[d]}
            assert image == stratum
            assert len(stratum) == len(primitives[d])
        orbits = {(l, min(rotate(ws, d) for d in range(n)))
                  for l, ws in primitives[n]}
        assert len(primitives[n]) == n * len(orbits)
        assert len(orbits) == ledger[n]
    print(f"word tuples: r={r}, 1<=n<={bound}; counts, root bijections, actual orbit sizes OK")


def main():
    bound = 20
    published = {
        2: [1, 1, 3, 18, 213, 4128, 122638, 5096305, 284192429, 20375905738,
            1829560187405, 200829815300994, 26471873341135571,
            4124649654997542447, 750006492020987263020, 157382918361825037892997],
        3: [1, 1, 5, 86, 4052, 400401, 71827456, 21068995258, 9429303819612,
            6105894632883407, 5493030296624140330, 6644655430011095138676,
            10523095865317003368417750, 21337870239129956669159151372],
    }
    published_b = {2: [1, 1, 5, 46, 773, 19426, 708062, 34740805],
                   3: [1, 1, 9, 244, 15833, 1980126, 428447592, 146966837193]}
    ledgers = {}
    for r in range(1, 6):
        b = power_sums(r, bound)
        assert weighted_partition_product(r, bound) == b
        a = source_recurrence(b)
        assert a == formal_exponential(b)
        c = primitive_ledger(b)
        assert a == integral_euler_product(c)
        assert all(v.denominator == 1 and v >= 0 for v in a)
        assert a[0] == 1
        if r in published:
            assert a[:len(published[r])] == published[r]
            assert b[:len(published_b[r])] == published_b[r]
            print(f"source fidelity: r={r}; a[0..{len(published[r])-1}] and b[0..7] match OEIS")
        print(f"r={r}: b[1..5]={b[1:6]}; c[1..5]={c[1:6]}; a[0..5]={[int(v) for v in a[:6]]}")
        print(f"r={r}: exact partition product, exponential, recurrence and integral Euler product agree through n={bound}")
        ledgers[r] = c
    for r, limit in [(1, 6), (2, 5), (3, 4)]:
        check_roots_and_orbits(r, limit, ledgers[r])
    print("Finite consistency only; all-index proof remains a separate Lean obligation.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
