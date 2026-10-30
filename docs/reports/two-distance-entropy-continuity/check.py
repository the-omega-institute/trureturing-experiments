#!/usr/bin/env python3
"""Exact finite checks for Jabbour–Datta, arXiv:2310.17329v3.

All probabilities use Fraction. A logarithmic expression is represented as
rational coefficients of log_2(prime). Equality uses unique factorization;
order uses comparison of integer products after clearing denominators.
No approximate logarithms, floating point, or numerical tolerances are used.
These checks do not prove the universal quantum conjecture.
"""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import lcm


@lru_cache(None)
def factor(n):
    assert isinstance(n, int) and n > 0
    result = []
    p = 2
    while p * p <= n:
        count = 0
        while n % p == 0:
            n //= p
            count += 1
        if count:
            result.append((p, count))
        p += 1
    if n > 1:
        result.append((n, 1))
    return tuple(result)


def add(*terms):
    out = defaultdict(Q)
    for term in terms:
        for prime, coefficient in term.items():
            out[prime] += coefficient
    return {p: c for p, c in out.items() if c}


def scale(a, term):
    return {p: a * c for p, c in term.items() if a * c}


def log2(x):
    x = Q(x)
    assert x > 0
    return add(dict(factor(x.numerator)),
               {p: -Q(c) for p, c in factor(x.denominator)})


def xlogx(x):
    x = Q(x)
    assert x >= 0
    return {} if x == 0 else scale(x, log2(x))


def sign(term):
    denominator = lcm(*(c.denominator for c in term.values()))
    positive = negative = 1
    for p, c in term.items():
        exponent = int(c * denominator)
        if exponent >= 0:
            positive *= p ** exponent
        else:
            negative *= p ** (-exponent)
    return (positive > negative) - (positive < negative)


def entropy(p):
    assert all(x >= 0 for x in p) and sum(p) == 1
    return scale(-1, add(*(xlogx(x) for x in p)))


def binary_entropy(e):
    return entropy((e, 1 - e))


def floor(x):
    return x.numerator // x.denominator


def ceil(x):
    return -((-x.numerator) // x.denominator)


def parameters(d, e, v):
    assert 0 < v <= e <= 1
    k = floor(e / v)
    mu = e - k * v
    m = ceil(e / v)
    assert 0 <= mu < v and m == k + (mu > 0)
    assert 2 * m <= d
    return k, mu, m, d - m


def bound(d, e, v):
    if e == 0:
        assert v == 0
        return {}
    k, mu, _, r = parameters(d, e, v)
    return add(binary_entropy(e), scale(k, xlogx(v)), xlogx(mu),
               scale(e, log2(r)), scale(-1, xlogx(e)))


def distances(p, q):
    difference = tuple(a - b for a, b in zip(p, q))
    return sum(map(abs, difference), Q(0)) / 2, max(map(abs, difference))


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, length - 1):
                yield (first,) + rest


def source_extremizers(d, e, v):
    k, mu, m, r = parameters(d, e, v)
    # Equations (42)–(44), with their separate zero-remainder branch.
    capped = (v / e,) * k + ((mu / e,) if mu else ())
    assert len(capped) == m and sum(capped) == 1
    q = capped + (Q(0),) * r
    p = tuple((1 - e) * x for x in capped) + (e / r,) * r
    return p, q


def fidelity():
    count = 0
    # Theorem 3, all unordered pairs in these entire rational simplexes.
    for d in range(2, 6):
        for n in range(1, 7):
            states = [tuple(Q(x, n) for x in row)
                      for row in compositions(n, d)]
            entropies = [entropy(p) for p in states]
            for i, j in combinations_with_replacement(range(len(states)), 2):
                p, q = states[i], states[j]
                e, v = distances(p, q)
                difference = add(entropies[i], scale(-1, entropies[j]))
                absolute_difference = scale(sign(difference), difference)
                assert sign(add(bound(d, e, v),
                                scale(-1, absolute_difference))) >= 0, (d, n, p, q)
                count += 1
    print('FIDELITY: Theorem 3, d=2..5, denominator n=1..6,')
    print('  every unordered simplex pair including identical states:', count)


def settlement_identities():
    count = integer_ratio = noninteger_ratio = endpoint = 0
    for d in range(2, 11):
        for n in range(1, 25):
            for a in range(1, n + 1):
                for b in range(1, a + 1):
                    e, v = Q(a, n), Q(b, n)
                    if 2 * ceil(e / v) > d:
                        continue
                    k, mu, m, r = parameters(d, e, v)
                    beta = v / e
                    # Core identity for applying Theorem 3 to orthogonal Jordan states.
                    assert add(binary_entropy(e), scale(e, bound(d, Q(1), beta))) == bound(d, e, v)
                    p, q = source_extremizers(d, e, v)
                    assert len(p) == len(q) == d
                    assert all(x >= 0 for x in p + q)
                    assert sum(p) == sum(q) == 1
                    assert distances(p, q) == (e, v)
                    assert add(entropy(p), scale(-1, entropy(q))) == bound(d, e, v)
                    assert sign(bound(d, e, v)) >= 0
                    delta = tuple(x - y for x, y in zip(p, q))
                    plus = tuple(max(x, Q(0)) / e for x in delta)
                    minus = tuple(max(-x, Q(0)) / e for x in delta)
                    assert sum(plus) == sum(minus) == 1
                    assert all(x * y == 0 for x, y in zip(plus, minus))
                    assert distances(plus, minus) == (1, beta)
                    # Integer-ratio reduction to the source's Sason expression.
                    if mu == 0:
                        coarse = add(binary_entropy(e), scale(e, log2(beta * d - 1)))
                        assert bound(d, e, v) == coarse
                        integer_ratio += 1
                    else:
                        noninteger_ratio += 1
                    endpoint += (e == 1)
                    count += 1
    assert integer_ratio and noninteger_ratio and endpoint
    assert bound(2, Q(0), Q(0)) == {}
    print('SETTLEMENT IDENTITY + SOURCE SHARPNESS: d=2..10,')
    print('  e=a/n, v=b/n, 1<=b<=a<=n<=24, 2*ceil(e/v)<=d:', count)
    print('  integer ratio:', integer_ratio, '; noninteger ratio:', noninteger_ratio)
    print('  e=1 endpoint:', endpoint, '; equal-state extension: checked')


def negative_controls():
    # A noninteger ratio where both the remainder term and support rounding matter.
    d, e, v = 6, Q(3, 5), Q(1, 4)
    k, mu, m, r = parameters(d, e, v)
    assert mu > 0 and m == k + 1
    correct = bound(d, e, v)
    omitted_remainder = add(correct, scale(-1, xlogx(mu)))
    wrong_support = add(correct, scale(-e, log2(r)), scale(e, log2(d - k)))
    assert correct != omitted_remainder
    assert correct != wrong_support
    assert sign(log2(Q(3, 2))) == 1
    assert sign(log2(Q(2, 3))) == -1
    print('NEGATIVE CONTROLS: omitted remainder and floor-for-ceiling both detected')


if __name__ == '__main__':
    fidelity()
    settlement_identities()
    negative_controls()
    print('ALL_CHECKS_OK')
