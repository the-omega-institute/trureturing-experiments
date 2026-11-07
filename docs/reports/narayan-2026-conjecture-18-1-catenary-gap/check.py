#!/usr/bin/env python3
"""Exact finite checks for Narayan, arXiv:2610.02495v1, Conjecture 18.1.

Enumerate all nonnegative coordinate factorizations, retaining precisely
(a, ah+d, ..., ah+kd). Compute each catenary degree by Kruskal's algorithm
on the complete factorization-distance graph. Independently verify threshold
connectivity at the computed degree and disconnection one unit below it.
No floating point and no conjectured degree formula enters the computation.
"""
from collections import Counter
from itertools import combinations
from math import gcd


def factorization_table(generators, bound):
    """All coordinate vectors of value <= bound, by recursive enumeration."""
    table = [[] for _ in range(bound + 1)]
    z = [0] * len(generators)

    def visit(i, value):
        if i == len(generators):
            table[value].append(tuple(z))
            return
        weight = generators[i]
        assert weight > 0
        for count in range((bound - value) // weight + 1):
            z[i] = count
            visit(i + 1, value + count * weight)

    visit(0, 0)
    return table


def distance(z, w):
    common = tuple(min(x, y) for x, y in zip(z, w))
    return max(sum(x - c for x, c in zip(z, common)),
               sum(y - c for y, c in zip(w, common)))


def connected(size, edges, threshold, vertices=None):
    vertices = set(range(size)) if vertices is None else set(vertices)
    if len(vertices) <= 1:
        return True
    adjacency = {i: [] for i in vertices}
    for weight, i, j in edges:
        if weight <= threshold and i in vertices and j in vertices:
            adjacency[i].append(j)
            adjacency[j].append(i)
    start = next(iter(vertices))
    seen, todo = {start}, [start]
    while todo:
        for j in adjacency[todo.pop()]:
            if j not in seen:
                seen.add(j)
                todo.append(j)
    return seen == vertices


def exact_catenary(factorizations):
    size = len(factorizations)
    assert size > 0
    edges = [(distance(z, w), i, j)
             for (i, z), (j, w) in combinations(enumerate(factorizations), 2)]
    parent = list(range(size))

    def root(i):
        while i != parent[i]:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    degree, components = 0, size
    for weight, i, j in sorted(edges):
        u, v = root(i), root(j)
        if u != v:
            parent[u] = v
            components -= 1
            degree = weight
            if components == 1:
                break
    assert connected(size, edges, degree)
    if degree:
        assert not connected(size, edges, degree - 1)
    else:
        assert size == 1
    return degree, edges


def euclidean_values(x, y):
    assert x > 0 and y > 0
    result = {x, y}
    while x != y:
        if x > y:
            x -= y
        else:
            y -= x
        result.update((x, y))
    return result


def fidelity():
    """Reproduce Section 4 arithmetic and a Section 17.1 catenary witness."""
    expected = {1, 2, 3, 4, 7, 10, 13, 23, 33, 43, 76, 109, 142, 251, 393}
    arithmetic_tail = {393 - 142 * r for r in range(3)}
    assert arithmetic_tail == {393, 251, 109}
    assert euclidean_values(109, 142) | arithmetic_tail == expected
    bound = 12000
    table = factorization_table((4, 681, 790), bound)
    observed = set()
    count = 0
    for factorizations in table:
        count += len(factorizations)
        lengths = sorted({sum(z) for z in factorizations})
        observed.update(y - x for x, y in zip(lengths, lengths[1:]))
    assert observed <= expected, (observed, expected)
    print('FIDELITY source_section=4 generators=(4,681,790) bound=' + str(bound))
    print('FIDELITY factorizations=' + str(count) + ' observed_delta=' + str(sorted(observed)))
    print('FIDELITY scope=exact_tail_and_Euclidean_union; finite_delta_subset_only')
    witness = factorization_table((5, 11, 12, 13), 23)[23]
    degree, _ = exact_catenary(witness)
    assert set(witness) == {(2, 0, 0, 1), (0, 1, 1, 0)}
    assert degree == 3
    print('FIDELITY source_theorem=17.1 parameters=(5,2,1,3) n=23 c(n)=3=h+1')


def main():
    fidelity()
    counts = Counter()
    histogram = Counter()
    for a in range(5, 13):
        for k in range(3, a - 1):
            for h in range(2, 6):
                for d in range(1, 9):
                    if gcd(a, d) != 1:
                        continue
                    generators = (a,) + tuple(a * h + i * d for i in range(1, k + 1))
                    n2 = 2 * a * h + (k + 1) * d
                    nh = 2 * a * h + 3 * d
                    bound = max(n2, nh) + a * h
                    table = factorization_table(generators, bound)
                    counts['parameter_tuples'] += 1
                    observed = set()
                    for n, factorizations in enumerate(table):
                        if not factorizations:
                            continue
                        assert len(factorizations) == len(set(factorizations))
                        assert all(sum(x * g for x, g in zip(z, generators)) == n
                                   for z in factorizations)
                        degree, edges = exact_catenary(factorizations)
                        observed.add(degree)
                        histogram[degree] += 1
                        counts['represented_values'] += 1
                        counts['factorizations'] += len(factorizations)
                        counts['distance_pairs'] += len(edges)
                        assert degree in (0, 2) or degree >= h + 1, (a, h, d, k, n, degree)
                        m = [sum(z[1:]) for z in factorizations]
                        t = [sum(i * z[i] for i in range(1, k + 1)) for z in factorizations]
                        m0 = min(m)
                        block = {i for i in range(len(m)) if m[i] == m0}
                        assert all(m0 <= t[i] < m0 + a for i in block)
                        assert len({(t[i], factorizations[i][0]) for i in block}) == 1
                        assert connected(len(m), edges, 2, block)
                        counts['minimum_count_blocks'] += 1
                        for weight, i, j in edges:
                            if (i in block) != (j in block):
                                assert weight >= h + 1
                                counts['crossing_pairs'] += 1
                    assert observed.intersection(range(h + 2)) == {0, 2, h + 1}
                    c2, _ = exact_catenary(table[n2])
                    ch, _ = exact_catenary(table[nh])
                    assert c2 == 2
                    assert ch == h + 1
                    assert len(table[nh]) == 2
                    counts['exact_witness_checks'] += 2
    print('WINDOW a=5..12 h=2..5 d=1..8 gcd(a,d)=1 3<=k<=a-2')
    print('WINDOW n=0..max(2*a*h+(k+1)*d,2*a*h+3*d)+a*h; represented_values_only')
    for key in sorted(counts):
        print(key + '=' + str(counts[key]))
    print('catenary_degree_histogram=' + str(sorted(histogram.items())))
    print('CHECK method=complete_graph_Kruskal_plus_independent_threshold_BFS')
    print('CHECK scope=finite_only_not_a_universal_proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
