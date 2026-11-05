#!/usr/bin/env python3
"""Exact finite checks for Kubo's assignment coordinates and an inverse bound.

Source: Susumu Kubo, arXiv:2606.30184v1, equation (1), Section 3,
and Section 4's parity-box example. No floating point is used.
This finite experiment does not prove a statement over arbitrary real matrices.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
from random import Random


@lru_cache(None)
def assignment_groups(n, r):
    """All row choices: -1 means unused; a nonnegative label is a column."""
    groups = defaultdict(list)
    for labels in product(range(-1, r), repeat=n):
        c = tuple(labels.count(a) for a in range(r))
        groups[c].append(tuple((i, a) for i, a in enumerate(labels) if a >= 0))
    expected = {c for c in product(range(n + 1), repeat=r) if sum(c) <= n}
    assert set(groups) == expected
    return dict(groups)


def coordinates(M):
    n, r = len(M), len(M[0])
    return {c: max(sum((M[i][a] for i, a in edges), Q(0)) for edges in choices)
            for c, choices in assignment_groups(n, r).items()}


def subset_definition(M, c):
    """Independent literal disjoint-subsets evaluation of source equation (1)."""
    n, r = len(M), len(M[0])
    subsets = [tuple(i for i in range(n) if mask & (1 << i)) for mask in range(1 << n)]
    values = []
    for sets in product(subsets, repeat=r):
        if tuple(map(len, sets)) != c:
            continue
        flat = tuple(i for s in sets for i in s)
        if len(set(flat)) == len(flat):
            values.append(sum((M[i][a] for a, s in enumerate(sets) for i in s), Q(0)))
    assert values, (M, c)
    return max(values)


def G(M, lam):
    return sum((max(Q(0), max(x - t for x, t in zip(row, lam))) for row in M), Q(0))


def orthant(M, lam):
    return sum(all(x <= t for x, t in zip(row, lam)) for row in M)


def distance_squared(M, N):
    return min(sum((M[i][a] - N[j][a]) ** 2
                   for i, j in enumerate(p) for a in range(len(M[0])))
               for p in permutations(range(len(M))))


def check_pair(M, N, bM, bN):
    n, r = len(M), len(M[0])
    diffs = [abs(bM[c] - bN[c]) for c in bM]
    e, l2sq = max(diffs), sum(x * x for x in diffs)
    d2 = distance_squared(M, N)
    factor_squared = 64 * n ** 3 * r
    assert d2 <= factor_squared * e * e, (M, N, e, d2)
    assert d2 <= factor_squared * l2sq
    assert l2sq <= sum(map(sum, bM)) * d2  # source forward Lipschitz bound
    if r == 1:
        assert d2 <= 4 * l2sq  # source's explicit r=1 inverse bound
    return e


def grid_check(M, N, bM, bN, eta):
    n, r = len(M), len(M[0])
    e = max(abs(bM[c] - bN[c]) for c in bM)
    assert eta > e >= 0
    thresholds, cluster_maps = [], []
    for a in range(r):
        values = sorted({row[a] for row in M + N})
        cuts, cmap, label, start = [], {}, 0, values[0]
        for k, x in enumerate(values):
            if k and x - values[k - 1] > 4 * eta:
                assert values[k - 1] - start < 8 * n * eta
                cuts.append(Q(x + values[k - 1]) / 2)
                start, label = x, label + 1
            cmap[x] = label
        assert values[-1] - start < 8 * n * eta
        thresholds.append([values[0] - 3 * eta] + cuts + [values[-1] + 3 * eta])
        cluster_maps.append(cmap)
    countM, countN = {}, {}
    for indices in product(*(range(len(t)) for t in thresholds)):
        lam = tuple(t[j] for t, j in zip(thresholds, indices))
        for row in M + N:
            assert all(abs(x - t) > 2 * eta for x, t in zip(row, lam))
        gm, gn = G(M, lam), G(N, lam)
        assert abs(gm - gn) <= e
        minus, plus = tuple(t - 2 * eta for t in lam), tuple(t + 2 * eta for t in lam)
        nm, nn = orthant(M, lam), orthant(N, lam)
        assert G(M, minus) - G(M, plus) == 4 * eta * (n - nm)
        assert G(N, minus) - G(N, plus) == 4 * eta * (n - nn)
        assert 4 * eta * abs(nm - nn) <= 2 * e
        assert nm == nn
        countM[indices], countN[indices] = nm, nn
    boxesM, boxesN = defaultdict(list), defaultdict(list)
    for cloud, boxes in [(M, boxesM), (N, boxesN)]:
        for i, row in enumerate(cloud):
            boxes[tuple(cmap[x] for cmap, x in zip(cluster_maps, row))].append(i)
    for cell in product(*(range(len(t) - 1) for t in thresholds)):
        # A high endpoint is selected by bit 1; a low endpoint contributes -1.
        corners = [(tuple(j + bit for j, bit in zip(cell, bits)),
                    (-1) ** (r - sum(bits))) for bits in product((0, 1), repeat=r)]
        assert sum(sign * countM[idx] for idx, sign in corners) == len(boxesM[cell])
        assert sum(sign * countN[idx] for idx, sign in corners) == len(boxesN[cell])
        assert len(boxesM[cell]) == len(boxesN[cell])
    pairing = []
    for cell, rows in boxesM.items():
        pairing.extend(zip(rows, boxesN[cell]))
    assert sorted(i for i, _ in pairing) == list(range(n))
    assert sorted(j for _, j in pairing) == list(range(n))
    for i, j in pairing:
        assert all(abs(x - y) < 8 * n * eta for x, y in zip(M[i], N[j]))
    d2 = sum((M[i][a] - N[j][a]) ** 2 for i, j in pairing for a in range(r))
    assert d2 < 64 * n ** 3 * r * eta * eta


def main():
    A = ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0))
    B = ((0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 1, 1))
    a, b = coordinates(A), coordinates(B)
    assert (a[(2, 1, 1)], b[(2, 1, 1)], a[(1, 1, 1)], b[(1, 1, 1)]) == (3, 4, 3, 3)
    for cols in [(0, 1), (0, 2), (1, 2)]:
        assert Counter(tuple(row[j] for j in cols) for row in A) == Counter(
            tuple(row[j] for j in cols) for row in B)
    for M, bc in [(A, a), (B, b)]:
        for c in bc:
            assert bc[c] == subset_definition(M, c)
    print('FIDELITY parity-box: b_(2,1,1)=(3,4); b_(1,1,1)=(3,3); pairwise projections equal')
    pair_checks, transform_checks, grid_checks = 0, 0, 0
    configs = [(n, r, (-1, 0, 1)) for n, r in
               [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (3, 1), (3, 2)]]
    configs += [(n, 3, (0, 1)) for n in (2, 3)]
    for n, r, alphabet in configs:
        clouds = list(combinations_with_replacement(list(product(alphabet, repeat=r)), n))
        bcs = [coordinates(M) for M in clouds]
        for M, bc in zip(clouds, bcs):
            assert bc[(0,) * r] == 0
            for lam in product((Q(-3, 2), Q(0), Q(3, 2)), repeat=r):
                dual = max(value - sum(ca * t for ca, t in zip(c, lam)) for c, value in bc.items())
                assert G(M, lam) == dual
                transform_checks += 1
        local_pairs = 0
        for i, M in enumerate(clouds):
            for j in range(i, len(clouds)):
                N = clouds[j]
                e = check_pair(M, N, bcs[i], bcs[j])
                grid_check(M, N, bcs[i], bcs[j], e + Q(1, 100))
                grid_checks += 1
                pair_checks += 1
                local_pairs += 1
        print(f'EXHAUSTIVE n={n} r={r} alphabet={alphabet}: {len(clouds)} multisets, {local_pairs} unordered pairs')
    rng = Random(260630184)
    for n in range(1, 5):
        for r in range(1, 4):
            for trial in range(25):
                M = tuple(tuple(Q(rng.choice((-100, 0, 100))) for _ in range(r)) for _ in range(n))
                N = tuple(tuple(x + Q(rng.randrange(-2, 3), 1000) for x in row) for row in M)
                if trial == 0:
                    N = tuple(reversed(M))  # exact zero-error multiplicities
                else:
                    N = tuple(reversed(N))
                bcM, bcN = coordinates(M), coordinates(N)
                e = check_pair(M, N, bcM, bcN)
                for eta in (e + Q(1, 100), e + 1, 2 * e + Q(1, 17)):
                    grid_check(M, N, bcM, bcN, eta)
                    grid_checks += 1
                pair_checks += 1
    # A cut occurs only for a gap strictly larger than 4 eta, including equality cases.
    M, N = ((Q(0),), (Q(4),)), ((Q(0),), (Q(4),))
    grid_check(M, N, coordinates(M), coordinates(N), Q(1))
    grid_checks += 1
    print(f'CHECK_COUNTS pairs={pair_checks} transform_values={transform_checks} gap_grids={grid_checks}')
    print('EXACT_ARITHMETIC integers and fractions.Fraction only; square-root inequalities checked after squaring nonnegative quantities')
    print('SCOPE finite checks support the planned general argument; no proof of the universal claim is asserted')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
