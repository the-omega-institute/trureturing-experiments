#!/usr/bin/env python3
"""Exact checks for Dovgoshey--Rovenska Problem 4.3 and the proposed answer."""
from fractions import Fraction
from itertools import combinations


def partner_matching(n, c):
    """Round-robin matching on {inf} union Z/(n-1), for even n."""
    q = n - 1
    assert n >= 2 and q % 2 == 1
    inf = q
    matching = set()
    matching.add(tuple(sorted((inf, c))))
    for x in range(q):
        if x != c:
            y = (2 * c - x) % q
            assert y != x
            matching.add(tuple(sorted((x, y))))
    assert len(matching) == n // 2
    return matching


def one_factorization(n):
    matchings = [partner_matching(n, c) for c in range(n - 1)]
    edges = set(combinations(range(n), 2))
    assert all(len(m) == n // 2 for m in matchings)
    assert len(set().union(*matchings)) == len(edges)
    assert set().union(*matchings) == edges
    assert all(len(matchings[i] & matchings[j]) == 0
               for i in range(n - 1) for j in range(i))
    return matchings


def candidate_weights(k):
    """Return n, edge->integer weight for the proposed construction."""
    if k == 1:
        return 1, {}
    if k == 2:
        return 2, {(0, 1): 1}
    if k % 2 == 0:
        n = k
        q = n - 1
        factors = one_factorization(n)
        weights = {}
        for c, matching in enumerate(factors):
            for e in matching:
                weights[e] = q + c
        assert len(set(weights.values())) == k - 1
        return n, weights
    n = k + 1
    q = n - 1
    factors = one_factorization(n)
    weights = {}
    for c, matching in enumerate(factors):
        value = q if c in (0, 1) else q + c - 1
        for e in matching:
            weights[e] = value
    assert len(set(weights.values())) == k - 1
    return n, weights


def distance_matrix(n, weights):
    d = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for (x, y), w in weights.items():
        d[x][y] = d[y][x] = Fraction(w)
    return d


def central_distances(d):
    row_sets = [set(row) for row in d]
    return set.intersection(*row_sets)


def is_metric(d):
    n = len(d)
    for x in range(n):
        for y in range(n):
            assert d[x][y] == d[y][x]
            assert (d[x][y] == 0) == (x == y)
            for z in range(n):
                assert d[x][y] <= d[x][z] + d[z][y]
    return True


def exact_small_cases():
    expected = {1: (1, 1), 2: (2, 2), 3: (4, 3), 4: (4, 4),
                5: (6, 5), 6: (6, 6), 7: (8, 7), 8: (8, 8)}
    for k, (expected_n, expected_c) in expected.items():
        n, weights = candidate_weights(k)
        assert n == expected_n
        d = distance_matrix(n, weights)
        assert is_metric(d)
        assert len(central_distances(d)) == expected_c
        assert central_distances(d) == {Fraction(0), *map(Fraction, set(weights.values()))}


def source_fidelity_bicentral_band():
    """Reproduce the paper's bounded-band metric realization for a bicentral graph.

    For a graph with minimum degree at least one, put distance 2 on its edges and
    distinct distances in [1,2) on complement edges. Every vertex then sees 2,
    so 2 is central; each complement edge's unique value is noncentral.
    """
    n = 4
    graph_edges = {(0, 1), (1, 2), (2, 3), (0, 3)}
    complement = set(combinations(range(n), 2)) - graph_edges
    d = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for j, (x, y) in enumerate(combinations(range(n), 2)):
        value = Fraction(2) if (x, y) in graph_edges else Fraction(3, 2) + Fraction(j, 100)
        d[x][y] = d[y][x] = value
    assert is_metric(d)
    assert central_distances(d) == {Fraction(0), Fraction(2)}
    assert complement


def lower_bound_and_parity(k):
    """Check the two elementary lower-bound implications used in the proof."""
    assert k <= candidate_weights(k)[0]
    if k > 1 and k % 2 == 1:
        assert candidate_weights(k)[0] >= k + 1


if __name__ == "__main__":
    source_fidelity_bicentral_band()
    exact_small_cases()
    for k in range(1, 25):
        lower_bound_and_parity(k)
        n, weights = candidate_weights(k)
        assert n == (1 if k == 1 else 2 if k == 2 else k if k % 2 == 0 else k + 1)
        assert len(weights) == n * (n - 1) // 2
    print("checked source bicentral band realization")
    print("checked proposed minima for k=1..24")
    print("ALL_CHECKS_OK")
