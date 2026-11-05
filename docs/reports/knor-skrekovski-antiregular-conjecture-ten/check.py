#!/usr/bin/env python3
"""Exact degree-vector checks; rational root enclosures use integer powers only."""

import hashlib
import json
from collections import Counter
from pathlib import Path

SCALE = 1 << 96


def root_interval(a, exponent):
    """Return L,U with L/SCALE < a**(1/exponent) < U/SCALE."""
    target = a * SCALE**exponent
    lo, hi = SCALE, 2 * SCALE
    assert lo**exponent < target < hi**exponent
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**exponent <= target:
            lo = mid
        else:
            hi = mid
    assert lo**exponent < target < hi**exponent
    return lo, hi


def coefficients(a):
    a1, a2, a3, a4 = a
    return a1*a2 + a2*a3 + a3*a4, a1*a3 + a2*a4, a1*a4


def patterns(n):
    k, r = divmod(n, 4)
    return {
        0: [(k, k, k, k)],
        1: [(k, k, k, k+1)],
        2: [(k, k+1, k, k+1), (k+1, k, k+1, k)],
        3: [(k+1, k, k+1, k+1)],
    }[r]


def lower_difference(p, q, b, c):
    d = [x-y for x, y in zip(p, q)]
    return d[0]*SCALE + d[1]*b[0 if d[1] >= 0 else 1] + d[2]*c[0 if d[2] >= 0 else 1]


def check_vectors(n, exponent, check_identity=False):
    b, c = root_interval(2, exponent), root_interval(3, exponent)
    wanted = patterns(n)
    best = coefficients(wanted[0])
    assert all(coefficients(p) == best for p in wanted)
    count = 0
    minimum_gap_numerator = None
    for a1 in range(n+1):
        for a2 in range(n-a1+1):
            for a3 in range(n-a1-a2+1):
                a4 = n-a1-a2-a3
                if (a1+a3) % 2:
                    continue
                a = a1, a2, a3, a4
                count += 1
                score = coefficients(a)
                if check_identity:
                    u, x, y = a1+a4, a1-a4, a2-a3
                    identity = (n*n-u*u-y*y+2*x*y,
                                2*n*u-2*u*u-2*x*y, u*u-x*x)
                    assert tuple(4*t for t in score) == identity
                    assert (x-u) % 2 == (y-(n-u)) % 2 == (x-y+n) % 4 == 0
                if a in wanted:
                    assert score == best
                else:
                    gap = lower_difference(best, score, b, c)
                    assert gap > 0, (n, exponent, a, gap)
                    minimum_gap_numerator = gap if minimum_gap_numerator is None else min(minimum_gap_numerator, gap)
    return count, minimum_gap_numerator


def havel_hakimi(a):
    remaining = [d for d, count in enumerate(a, 1) for _ in range(count)]
    edges = set()
    while any(remaining):
        order = sorted(range(len(remaining)), key=lambda v: (-remaining[v], v))
        v, degree = order[0], remaining[order[0]]
        remaining[v] = 0
        for w in order[1:degree+1]:
            assert remaining[w] > 0
            remaining[w] -= 1
            edges.add(tuple(sorted((v, w))))
    # A graphical sequence need not yield a connected realization under the
    # greedy algorithm. Join components by a degree-preserving two-switch.
    n = len(remaining)
    while True:
        parts = components(n, edges)
        if len(parts) == 1:
            return edges
        joined = False
        for A, B in sorted(edges):
            if len(components(n, edges - {(A, B)})) != len(parts):
                continue
            own = next(part for part in parts if A in part)
            C, D = next(e for e in sorted(edges) if e[0] not in own)
            edges = (edges - {(A, B), (C, D)}) | {
                tuple(sorted((A, C))), tuple(sorted((B, D)))}
            joined = True
            break
        assert joined


def components(n, edges):
    neighbors = [set() for _ in range(n)]
    for v, w in edges:
        neighbors[v].add(w)
        neighbors[w].add(v)
    unseen = set(range(n))
    parts = []
    while unseen:
        seen = {min(unseen)}
        queue = list(seen)
        for v in queue:
            for w in neighbors[v] - seen:
                seen.add(w)
                queue.append(w)
        unseen -= seen
        parts.append(seen)
    return parts


def graph_data(n, edges):
    assert len(edges) == len(set(edges))
    neighbors = [set() for _ in range(n)]
    for v, w in edges:
        assert 0 <= v < w < n
        neighbors[v].add(w)
        neighbors[w].add(v)
    visited = {0}
    queue = [0]
    for v in queue:
        for w in neighbors[v] - visited:
            visited.add(w)
            queue.append(w)
    assert len(visited) == n, (n, sorted(edges))
    degrees = [len(s) for s in neighbors]
    assert min(degrees) >= 1 and max(degrees) <= 4
    histogram = Counter(degrees)
    return degrees, tuple(histogram[d] for d in range(1, 5))


def insert(n, edges):
    A, B = min(edges)
    a, b, c, d = range(n, n+4)
    extra = [(A, d), (d, c), (c, B), (b, c), (b, d), (a, d)]
    return (edges - {(A, B)}) | {tuple(sorted(e)) for e in extra}


def main():
    # Theorem 9, n=8, f(8)=1/64. Its logarithmic hypothesis is
    # 3**(1/64) <= 24/23, certified here without logarithms.
    assert 3 * 23**64 < 24**64
    count, gap = check_vectors(8, 64)
    print(f'FIDELITY Theorem 9 n=8 exponent=1/64 vectors={count} strict_gap_lower={gap}/{SCALE}')
    print('ROOT_INTERVALS integer-power brackets; denominator=' + str(SCALE))
    total = 0
    for n in range(7, 65):
        count, gap = check_vectors(n, n, check_identity=True)
        total += count
        print(f'CLAIM n={n} vectors={count} maximizers={json.dumps(patterns(n))} strict_gap_lower={gap}/{SCALE}')
    print(f'IDENTITY_AND_CLAIM vectors_total={total}; all even-degree-sum vectors, including nonrealizable ones')
    bases = []
    for n in range(7, 11):
        for target in patterns(n):
            edges = havel_hakimi(target)
            degrees, histogram = graph_data(n, edges)
            assert histogram == target
            # Direct pair summation and degree-class aggregation agree.
            direct = [0, 0, 0]
            for v in range(n):
                for w in range(v+1, n):
                    difference = abs(degrees[v]-degrees[w])
                    if difference:
                        direct[difference-1] += 1
            assert tuple(direct) == coefficients(target)
            bases.append({'n': n, 'counts': target, 'edges': sorted(edges)})
            for step in range(25):
                size = n+4*step
                before, hist = graph_data(size, edges)
                edges = insert(size, edges)
                after, new_hist = graph_data(size+4, edges)
                assert after[:size] == before
                assert after[size:] == [1, 2, 3, 4]
                assert new_hist == tuple(x+1 for x in hist)
                assert new_hist in patterns(size+4)
    print('CONNECTED_BASES ' + json.dumps(bases, separators=(',', ':')))
    print('INSERTION_CHECKS bases=5 iterations_per_base=25 old_degrees_preserved new_degrees=[1,2,3,4] connected=true')
    print('CHECK_SHA256 ' + hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
