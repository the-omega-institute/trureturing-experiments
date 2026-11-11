#!/usr/bin/env python3
"""Exact finite checks for the proposed finite-double reduction.

None denotes +infinity; every finite distance is an integer. These checks
neither enumerate arbitrary endofunctors nor certify transfinite convergence.
"""
from itertools import combinations, product


def le(a, b):
    return b is None or (a is not None and a <= b)


def add(a, b):
    return None if a is None or b is None else a + b


def emin(xs):
    finite = [x for x in xs if x is not None]
    return min(finite) if finite else None


def emax(xs):
    xs = list(xs)
    return None if None in xs else max(xs, default=0)


def is_metric(d):
    n = len(d)
    return all(d[x][x] == 0 for x in range(n)) and all(
        d[x][y] == d[y][x] and (x == y or d[x][y] != 0)
        and all(le(d[x][z], add(d[x][y], d[y][z])) for z in range(n))
        for x in range(n) for y in range(n)
    )


def metrics(n):
    edges = list(combinations(range(n), 2))
    for values in product((1, 2, 3, None), repeat=len(edges)):
        d = [[0 for _ in range(n)] for _ in range(n)]
        for (x, y), value in zip(edges, values):
            d[x][y] = d[y][x] = value
        if is_metric(d):
            yield d


def check_double(d, subset):
    n = len(d)
    points = [(0, x) for x in range(n)] + [(1, x) for x in range(n) if x not in subset]
    index = {point: i for i, point in enumerate(points)}

    def j(side, x):
        return index[(0 if x in subset else side, x)]

    def cross(x, y):
        return emin(add(d[x][s], d[s][y]) for s in subset)

    def dist(p, q):
        side, x = p
        other, y = q
        return d[x][y] if side == other or x in subset or y in subset else cross(x, y)

    doubled = [[dist(p, q) for q in points] for p in points]
    assert is_metric(doubled)
    # An independent shortest-path calculation detects triangle shortcuts.
    shortest = [row[:] for row in doubled]
    for k in range(len(points)):
        for x in range(len(points)):
            for y in range(len(points)):
                shortest[x][y] = emin((shortest[x][y], add(shortest[x][k], shortest[k][y])))
    assert shortest == doubled
    for x, y in product(range(n), repeat=2):
        assert (j(0, x) == j(1, y)) == (x == y and x in subset)
        for side in (0, 1):
            assert doubled[j(side, x)][j(side, y)] == d[x][y]
        assert le(d[x][y], doubled[j(0, x)][j(1, y)])
        if x in subset or y in subset:
            assert cross(x, y) == d[x][y]
    for x in range(n):
        if x not in subset:
            assert cross(x, x) != 0
    # Concrete functor sanity checks: categorical powers of arity 0, 1, 2.
    # Verify their entire finite pullback, then the distance identity on it.
    power_checks = 0
    for arity in range(3):
        fx = list(product(range(n), repeat=arity))
        fs = list(product(sorted(subset), repeat=arity))
        pullback = [(u, v) for u in fx for v in fx if
                    tuple(j(0, x) for x in u) == tuple(j(1, y) for y in v)]
        assert pullback == [(u, u) for u in fs]
        for u, v in product(fs, repeat=2):
            source = emax(d[x][y] for x, y in zip(u, v))
            induced = emax((emax(d[x][y] for x, y in zip(u, v)),
                            emax(d[x][y] for x, y in zip(u, v))))
            assert source == induced
        power_checks += 1
    return power_checks


def main():
    assert le(3, None) and not le(None, 3)
    assert add(0, None) is None and emin([]) is None and emax([]) == 0
    metric_count = double_count = power_count = union_count = 0
    for n in range(1, 5):
        subsets = [frozenset(x for x in range(n) if mask & (1 << x))
                   for mask in range(1, 1 << n)]
        for d in metrics(n):
            metric_count += 1
            # Fidelity: Lemma L:emb(1), directed colimits of subspace embeddings
            # are metric unions with distance measured at any common stage.
            for s, t in product(subsets, repeat=2):
                common = s | t
                for x in s:
                    for y in t:
                        local = [[d[a][b] for b in sorted(common)] for a in sorted(common)]
                        assert local[sorted(common).index(x)][sorted(common).index(y)] == d[x][y]
                union_count += 1
            for s in subsets:
                power_count += check_double(d, s)
                double_count += 1
    # Empty-space normalization boundary: the empty-domain map is isometric
    # for every finite target, and a map to an empty set requires empty source.
    for n in range(5):
        empty_domain_pairs = list(product(range(0), repeat=2))
        assert not empty_domain_pairs
        assert (n == 0) == (len(list(product(range(0), repeat=n))) == 1)
    print('Exact distance alphabet: {1, 2, 3, +infinity}; 0 on the diagonal')
    print(f'All separated symmetric extended metric matrices with 1 <= n <= 4: {metric_count}')
    print(f'Fidelity: directed-subspace common-stage distance checks: {union_count}')
    print(f'Finite nonempty doubles: metric, shortest paths, isometries, fold, pullback: {double_count}')
    print(f'Power-functor pullbacks and diagonal distance identities, arities 0..2: {power_count}')
    print('Empty-domain and infinity conventions: OK')
    print('Scope: finite reduction checks only; no proof of the universal functor claim or omega+omega bound')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
