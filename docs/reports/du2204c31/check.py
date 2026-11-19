"""Exact checks for the structural coindex assertion in Du et al., Conjecture 3.1."""

from collections import Counter
from itertools import combinations
from math import isqrt


def initial_edge():
    return [{1}, {0}]


def attach(graph, a, b):
    """Adjoin one vertex to the endpoints of an existing edge."""
    assert b in graph[a] and a in graph[b]
    v = len(graph)
    graph.append({a, b})
    graph[a].add(v)
    graph[b].add(v)


def path_square(n):
    graph = initial_edge()
    for v in range(2, n):
        attach(graph, v - 2, v - 1)
    return graph


def book(n):
    graph = initial_edge()
    for _ in range(n - 2):
        attach(graph, 0, 1)
    return graph


def family(k):
    graph = path_square(2 * k)
    for i in range(k):
        for _ in range(k - 2):
            attach(graph, 2 * i, 2 * i + 1)
    return graph


def coindex_radicals(graph):
    """Count weights by integer radicand, once per unordered nonedge."""
    squares = [len(neighbors) ** 2 for neighbors in graph]
    return Counter(
        squares[u] + squares[v]
        for u, v in combinations(range(len(graph)), 2)
        if v not in graph[u]
    )


def scaled_interval(radicals, scale):
    lower = upper = 0
    for radicand, multiplicity in radicals.items():
        scaled_square = radicand * scale**2
        root_floor = isqrt(scaled_square)
        lower += multiplicity * root_floor
        upper += multiplicity * (root_floor + (root_floor**2 != scaled_square))
    return lower, upper


def check_fidelity():
    for n in range(5, 17):
        assert coindex_radicals(book(n)) == {8: (n - 2) * (n - 3) // 2}
        graph = book(n - 1)
        attach(graph, 0, 2)
        # In L_n, the three nonedge types need not have distinct radicands.
        expected = Counter()
        for radicand, multiplicity in [
            (8, (n - 3) * (n - 4) // 2),
            ((n - 2) ** 2 + 4, 1),
            (13, n - 4),
        ]:
            expected[radicand] += multiplicity
        assert coindex_radicals(graph) == expected
    print("Fidelity: Proposition 2.1 coindices of X_n and L_n, 5 <= n <= 16")


def check_uniform_inequalities():
    # Substitute k = 32 + t. The coefficients certify both inequalities for t >= 0.
    # 10(L-1)-9n = k^2-20k-10; 10(k-1)(k-2)-9n = k^2-30k+20.
    for linear, constant in [(-20, -10), (-30, 20)]:
        coefficients = (1, 64 + linear, 32**2 + 32 * linear + constant)
        assert all(coefficient >= 0 for coefficient in coefficients)
    # sqrt(2) < 3/2 and (2 - 81/100)*(3/2) < 9/5.
    assert 2 * 2**2 < 3**2
    assert (200 - 81) * 3 < 2 * 180
    print("Uniform algebra: both count bounds and a strict radical comparison for k >= 32")


def check_family(k):
    graph = family(k)
    n = k**2
    core = 2 * k
    private = n - core
    assert len(graph) == n
    assert sum(map(len, graph)) == 2 * (2 * n - 3)
    for v, neighbors in enumerate(graph):
        if v < core:
            assert len(neighbors) >= k
            assert private - sum(w >= core for w in neighbors) == (k - 1) * (k - 2)
        else:
            assert len(neighbors) == 2 and all(w < core for w in neighbors)
    figure = path_square(n)
    for v, neighbors in enumerate(figure):
        assert neighbors == set(range(max(0, v - 2), min(n, v + 3))) - {v}
    scale = 10**8
    witness_lower, _ = scaled_interval(coindex_radicals(graph), scale)
    _, figure_upper = scaled_interval(coindex_radicals(figure), scale)
    assert witness_lower > figure_upper
    print(f"k={k}, n={n}: coindex gap >= {witness_lower - figure_upper}/{scale} > 0")


def main():
    check_fidelity()
    check_uniform_inequalities()
    for k in (32, 33, 36, 40, 48):
        check_family(k)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
