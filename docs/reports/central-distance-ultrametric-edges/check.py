#!/usr/bin/env python3
"""Exact checks for arXiv:2610.11596v1, Theorems 3.10/3.13 and Conjecture 4.2.

Only integer and Fraction arithmetic is used. Finite checks are experimental
support for the stated unbounded argument, not a Lean proof.
"""
from fractions import Fraction as Q
from itertools import combinations, product


def require(test, message):
    if not test:
        raise AssertionError(message)


def pairs(n):
    return list(combinations(range(n), 2))


def table(n, labels):
    d = [[Q(0) for _ in range(n)] for _ in range(n)]
    for (i, j), v in labels.items():
        require(i < j, "unordered-pair order")
        d[i][j] = d[j][i] = Q(v)
    require(all(d[i][j] > 0 for i, j in pairs(n)), "missing positive distance")
    return d


def metric(d):
    n = len(d)
    return (
        all(d[i][j] == d[j][i] and d[i][j] >= 0
            and ((d[i][j] == 0) == (i == j))
            for i in range(n) for j in range(n))
        and all(d[i][j] <= d[i][k] + d[k][j]
                for i in range(n) for j in range(n) for k in range(n))
    )


def ultra(d):
    return metric(d) and all(
        d[i][j] <= max(d[i][k], d[k][j])
        for i in range(len(d)) for j in range(len(d)) for k in range(len(d))
    )


def center(d):
    return set.intersection(*(set(row) for row in d))


def central_edges(d):
    c = center(d)
    return {e for e in pairs(len(d)) if d[e[0]][e[1]] in c}


def distance_classes(d):
    return {r: {e for e in pairs(len(d)) if d[e[0]][e[1]] == r}
            for r in center(d) if r > 0}


def covers(n, edges):
    return set().union(*(set(e) for e in edges)) == set(range(n))


def star(vertices, root):
    return {tuple(sorted((root, v))) for v in vertices if v != root}


def is_star(n, edges):
    return any(edges == star(range(n), root) for root in range(n))


def condition_i_witness(d):
    """Check the specific partition A={0,1}, B={2,...,n-1}, star centers 0,2."""
    n = len(d)
    a, b = {0, 1}, set(range(2, n))
    diameter = max(max(row) for row in d)
    diam_edges = {e for e in pairs(n) if d[e[0]][e[1]] == diameter}
    cross = {(i, j) for i in a for j in b}
    return (min(len(a), len(b)) == 2 and diam_edges == cross
            and central_edges(d) == cross | star(a, 0) | star(b, 2))


def check_lower_bound(d):
    n = len(d)
    classes = distance_classes(d)
    require(len(center(d)) == 3 and len(classes) == 2, "three-central")
    for edges in classes.values():
        require(covers(n, edges), "central class is a spanning edge cover")
        require(2 * len(edges) >= n, "edge-cover incidence lower bound")
        require(len(edges) >= (n + 1) // 2, "edge-cover ceiling lower bound")
    a, b = classes.values()
    require(a.isdisjoint(b), "different distance classes disjoint")
    require(len(central_edges(d)) >= 2 * ((n + 1) // 2), "three-central bound")


def ultrametric_family(n):
    def distance(i, j):
        if i < 2 <= j:
            return 2
        if (i, j) == (0, 1) or i == 2:
            return 1
        return Q(1, 2)
    return table(n, {e: distance(*e) for e in pairs(n)})


def sparse_comparator(n):
    cycle = [tuple(sorted((i, (i + 1) % n))) for i in range(n)]
    a, b = set(cycle[::2]), set(cycle[1::2])
    if n % 2:
        require((0, 2) not in cycle, "odd chord is not a cycle edge")
        b.add((0, 2))
    require(a.isdisjoint(b) and covers(n, a) and covers(n, b), "two disjoint covers")
    require(len(a) == len(b) == (n + 1) // 2, "optimal cover cardinalities")
    remaining = [e for e in pairs(n) if e not in a | b]
    labels = {e: Q(4) for e in a} | {e: Q(5) for e in b}
    labels.update({e: Q(6) + Q(j, len(remaining) + 1)
                   for j, e in enumerate(remaining, 1)})
    return table(n, labels), a, b


def check_fidelity():
    """Exhaust all tables with positive distances in {1,2,3}, n=2,3,4."""
    for n in range(2, 5):
        es = pairs(n)
        metric_count = ultra_count = nonempty_count = 0
        minimum = n * n
        for values in product(range(1, 4), repeat=len(es)):
            d = table(n, dict(zip(es, values)))
            if not metric(d):
                continue
            metric_count += 1
            edges = central_edges(d)
            if edges:
                nonempty_count += 1
                minimum = min(minimum, len(edges))
                require(len(edges) >= (n + 1) // 2, "source Theorem 3.10")
            if len(center(d)) == 3:
                check_lower_bound(d)
            if ultra(d):
                ultra_count += 1
                require(len(edges) >= n - 1, "source Theorem 3.13 lower bound")
                require((len(edges) == n - 1) == is_star(n, edges),
                        "source Theorem 3.13 equality iff star")
        require(minimum == (n + 1) // 2, "source Theorem 3.10 sharpness")
        print(f"FIDELITY n={n}: metric_tables={metric_count}, "
              f"nonempty_central={nonempty_count}, ultrametric_tables={ultra_count}, "
              f"minimum_nonempty_central_edges={minimum}; Theorems 3.10,3.13 OK")


def check_four_points():
    x = [[Q(v) for v in row] for row in
         [[0, 1, 2, 2], [1, 0, 2, 2], [2, 2, 0, 1], [2, 2, 1, 0]]]
    y = [[Q(v) for v in row] for row in
         [[0, 4, 6, 5], [4, 0, 5, 7], [6, 5, 0, 4], [5, 7, 4, 0]]]
    require(ultra(x) and condition_i_witness(x), "four-point X, condition (i)")
    require(metric(y) and not ultra(y), "four-point Y is a non-ultrametric metric")
    require(center(x) == {0, 1, 2} and center(y) == {0, 4, 5}, "four-point centers")
    require(central_edges(x) == set(pairs(4)), "four-point X edges")
    expected = {(0, 1), (1, 2), (2, 3), (0, 3)}
    require(central_edges(y) == expected, "four-point Y edges")
    check_lower_bound(y)
    require(len(central_edges(x)) > len(central_edges(y)), "condition (ii) fails")
    print("FOUR_POINT X: center={0,1,2}, central_edges=6, condition_i=True")
    print("FOUR_POINT Y: center={0,4,5}, central_edges=4, metric=True, ultrametric=False")
    print("FOUR_POINT exact central Y edges:", sorted(central_edges(y)))


def check_families():
    for n in range(4, 33):
        x = ultrametric_family(n)
        require(ultra(x), "family X ultrametric")
        require(center(x) == {0, 1, 2}, "family X center")
        require(condition_i_witness(x), "family X condition (i)")
        require(len(central_edges(x)) == 3 * n - 6, "family X edge count")
        y, a, b = sparse_comparator(n)
        require(metric(y), "family comparator metric")
        require(not ultra(y), "family comparator is not ultrametric")
        require(center(y) == {0, 4, 5}, "family comparator exact center")
        require(central_edges(y) == a | b, "family comparator exact edges")
        require(len(central_edges(y)) == 2 * ((n + 1) // 2), "sharp comparator size")
        require(all(4 <= y[i][j] < 7 for i, j in pairs(n)), "metric band")
        check_lower_bound(y)
        require(3 * n - 6 > 2 * ((n + 1) // 2), "strict refuting gap")
        print(f"FAMILY n={n}: ultrametric_edges={3*n-6}, "
              f"metric_edges={2*((n+1)//2)}, centers=3, condition_i=True, "
              f"condition_ii=False, cover_sizes={len(a)},{len(b)}")
    print("SCOPE: finite family checks n=4..32; unbounded claims require the general proof.")


def main():
    check_fidelity()
    check_four_points()
    check_families()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
