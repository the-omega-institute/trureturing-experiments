#!/usr/bin/env python3
"""Exact checks for arXiv:2606.00420v3, Remark 4.11.

Independently written from the mathematical definitions; standard library only.
SPDX-License-Identifier: Apache-2.0
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import isqrt, lcm, prod


def square(x):
    x = Q(x)
    return x >= 0 and isqrt(x.numerator) ** 2 == x.numerator and isqrt(x.denominator) ** 2 == x.denominator


def square_class(n):
    assert n > 0
    result, prime = 1, 2
    while prime * prime <= n:
        parity = 0
        while n % prime == 0:
            n //= prime
            parity ^= 1
        if parity:
            result *= prime
        prime += 1
    return result * n


def factor(p, q):
    assert -1 < p < 1 and -1 < q < 1
    assert 1 + p*p > 0 and 1 - q*q > 0
    return (1 - p*p) / (1 + p*p) * (1 + q*q) / (1 - q*q)


def edge(p, q):
    a = factor(p, q)
    assert 1 + p*q > 0 and 1 - p*q > 0
    return a * (1 - p*q) / (1 + p*q), a * (1 + p*q) / (1 - p*q)


def normalize(parts):
    parts = tuple(map(Q, parts))
    assert all(x > 0 for x in parts)
    z = lcm(*(a for x in parts for a in (x.numerator, x.denominator)))
    result = tuple(z * x for x in parts)
    assert z > 0 and all(x.denominator == 1 for x in result)
    return z, result


def products(parts, k):
    return Counter(prod(c) for c in combinations(parts, k))


def pairings(indices):
    if not indices:
        yield ()
    else:
        first, *rest = indices
        for k, other in enumerate(rest):
            for tail in pairings(rest[:k] + rest[k+1:]):
                yield ((first, other),) + tail


def main():
    # Fidelity: Proposition 2.6's printed normalization example.
    source = tuple(map(Q, ('5/2', '12/5', '4/3', '4/7', '1/2', '7/16')))
    z, normalized = normalize(source)
    assert prod(source) == 1
    assert z == 1680
    assert normalized == (4200, 4032, 2240, 960, 840, 735)
    assert tuple(sorted((Q(z*z) / x for x in normalized), reverse=True)) == (3840, 3360, 2940, 1260, 700, 672)
    print('FIDELITY Proposition 2.6: Z=1680; Psi=(4200,4032,2240,960,840,735); flip_3 verified')

    # Fidelity: the first edge calculation in Example 3.6.
    assert edge(Q(3, 4), Q(1, 2)) == (Q(7, 33), Q(77, 75))
    print('FIDELITY Example 3.6: h+(3/4,1/2)=7/33; h-=77/75')

    # Fidelity: Remark 4.11's equal-size elementary-symmetric collision.
    target = (90, 50, 45, 25, 12, 12)
    partner = (75, 75, 36, 20, 18, 10)
    assert sum(target) == sum(partner) == 234
    assert products(target, 3) == products(partner, 3)
    assert prod(target) == 900**3
    assert tuple(sorted((Q(900, x) for x in target), reverse=True)) == partner
    print('FIDELITY Remark 4.11: sizes=234; all 20 triple-product occurrences agree; flip_3 verified')

    grid = sorted({Q(a, b) for b in range(1, 11) for a in range(-b+1, b)})
    scales = (1, 2, 3, 5, 7, 30, 1680)
    edge_count = 0
    for p in grid:
        for q in grid:
            hp, hm = edge(p, q)
            a = factor(p, q)
            assert hp > 0 and hm > 0 and hp * hm == a*a
            for scale in scales:
                assert (scale*hp) * (scale*hm) == (scale*a)**2
            edge_count += 1
    print(f'EDGE_IDENTITY: {len(grid)} rational coordinates; {edge_count} ordered pairs; {len(scales)} scales each; all exact')

    # All directed graphs on each prefix of a fixed four-vertex grid.
    vertices = (Q(-1, 2), Q(-1, 3), Q(0), Q(1, 2))
    tested, eulerian = 0, 0
    for n in range(5):
        candidates = [(u, v) for u in vertices[:n] for v in vertices[:n] if u != v]
        for mask in range(1 << len(candidates)):
            tested += 1
            edges = [e for i, e in enumerate(candidates) if mask & (1 << i)]
            if any(sum(u == v for u, _ in edges) != sum(w == v for _, w in edges) for v in vertices[:n]):
                continue
            eulerian += 1
            pairs = [edge(*e) for e in edges]
            parts = tuple(x for pair in pairs for x in pair)
            assert prod(parts) == 1
            scale, output = normalize(parts)
            assert Counter(output) != Counter(target)
            assert len(output) == 2*len(edges)
            for (p, q), (hp, hm) in zip(edges, pairs):
                assert (scale*hp)*(scale*hm) == (scale*factor(p, q))**2
    print(f'GRAPH_CHECK: {tested} loopless simple directed graphs; {eulerian} Eulerian graphs; all occurrence pairs and normalization verified')

    all_pairings = list(pairings(list(range(len(target)))))
    valid = [pairing for pairing in all_pairings if all(square(target[i]*target[j]) for i, j in pairing)]
    assert len(all_pairings) == 15 and not valid
    assert [i for i, x in enumerate(target) if square(x)] == [3]
    assert not any(square(target[3]*x) for i, x in enumerate(target) if i != 3)
    classes = tuple(square_class(x) for x in target)
    assert classes == (10, 2, 5, 1, 3, 3)
    eligible = [(i, j) for i, j in combinations(range(6), 2) if square(target[i]*target[j])]
    assert eligible == [(4, 5)]
    print('TARGET: square classes=(10,2,5,1,3,3); only compatible occurrence pair=(4,5), the two 12s')
    print('TARGET: 15 occurrence-level perfect pairings checked; valid=0; unique square occurrence is 25')

    rational_scales = sorted({Q(a, b) for a in range(1, 11) for b in range(1, 11)})
    for scale in rational_scales:
        assert not any(all(square(scale*target[i]*scale*target[j]) for i, j in pairing) for pairing in all_pairings)
    for multiple in range(1, 201):
        scaled_classes = [square_class(multiple*x) for x in target]
        assert any(count % 2 for count in Counter(scaled_classes).values())
    print(f'SCALING_CHECK: {len(rational_scales)} positive rational scales and integral multiples 1..200 retain the obstruction')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
