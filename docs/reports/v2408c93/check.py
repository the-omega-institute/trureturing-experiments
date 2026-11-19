#!/usr/bin/env python3
"""Check source examples and finite binary-incidence graphs with exact arithmetic."""

from fractions import Fraction
from itertools import combinations, permutations, product


def graph(order, edges):
    neighbors = [set() for _ in range(order)]
    for u, v in edges:
        assert u != v
        neighbors[u].add(v)
        neighbors[v].add(u)
    return neighbors


def complement(neighbors):
    vertices = set(range(len(neighbors)))
    return [vertices - row - {v} for v, row in enumerate(neighbors)]


def is_induced_path(neighbors, path):
    return len(set(path)) == len(path) and all(
        (path[j] in neighbors[path[i]]) == (j == i + 1)
        for i, j in combinations(range(len(path)), 2)
    )


def induced_path_supports(neighbors):
    def extend(path):
        yield frozenset(path)
        for vertex in neighbors[path[-1]]:
            if vertex not in path and not neighbors[vertex].intersection(path[:-1]):
                yield from extend((*path, vertex))

    return set().union(*(set(extend((v,))) for v in range(len(neighbors))))


def is_independent(neighbors, vertices):
    return all(v not in neighbors[u] for u, v in combinations(vertices, 2))


def by_definition(neighbors):
    """Enumerate color functions and test their fibers against all induced paths."""
    vertices = range(len(neighbors))
    supports = induced_path_supports(neighbors)
    for count in range(len(neighbors) + 1):
        for coloring in product(range(count), repeat=len(neighbors)):
            classes = [{v for v in vertices if coloring[v] == c} for c in range(count)]
            if all(
                is_independent(neighbors, cls)
                and all(len(cls & path) <= 2 for path in supports)
                for cls in classes
            ):
                return count
    raise AssertionError("singleton coloring missing")


def connected(neighbors):
    reached = {0}
    frontier = [0]
    while frontier:
        for v in neighbors[frontier.pop()] - reached:
            reached.add(v)
            frontier.append(v)
    return len(reached) == len(neighbors)


def incidence_graph(k):
    words = list(product((0, 1), repeat=k))
    tests = list(product(range(k), (0, 1)))
    r = len(words)
    return graph(
        r + len(tests),
        (
            (x, r + t)
            for x, word in enumerate(words)
            for t, (i, bit) in enumerate(tests)
            if word[i] == bit
        ),
    ), words, tests


def triple_path(neighbors, triple):
    """Find a P5 in the adjacency lists whose word vertices are this triple."""
    for left, middle, right in permutations(triple):
        first_tests = (neighbors[left] & neighbors[middle]) - neighbors[right]
        second_tests = (neighbors[middle] & neighbors[right]) - neighbors[left]
        if first_tests and second_tests:
            return left, min(first_tests), middle, min(second_tests), right
    raise AssertionError(f"no induced P5 for {triple}")


def verify_pair_coloring(neighbors, classes):
    vertices = [v for cls in classes for v in cls]
    assert sorted(vertices) == list(range(len(neighbors)))
    assert all(1 <= len(cls) <= 2 and is_independent(neighbors, cls) for cls in classes)


def check_incidence(k):
    neighbors, words, tests = incidence_graph(k)
    other = complement(neighbors)
    r = len(words)
    assert r == 2**k and len(tests) == 2 * k
    assert len(neighbors) == r + 2 * k
    assert connected(neighbors) and connected(other)

    checked = 0
    for triple in combinations(range(r), 3):
        path = triple_path(neighbors, triple)
        assert set(path[::2]) == set(triple) and is_induced_path(neighbors, path)
        checked += 1
    assert all(v in other[u] for u, v in combinations(range(r), 2))

    word_id = {word: i for i, word in enumerate(words)}
    half = r // 2
    original_classes = [(i, i + half) for i in range(half)]
    original_classes += [(r + 2 * i, r + 2 * i + 1) for i in range(k)]
    verify_pair_coloring(neighbors, original_classes)

    other_classes = [[i] for i in range(r)]
    assigned = []
    for t, (coord, bit) in enumerate(tests):
        unit = tuple(int(i == coord) for i in range(k))
        word = unit if bit else tuple(1 - b for b in unit)
        assigned.append(word_id[word])
        other_classes[word_id[word]].append(r + t)
    assert len(set(assigned)) == 2 * k
    verify_pair_coloring(other, other_classes)
    assert len(other_classes) == r and len(original_classes) == half + k

    lower = half * r
    conjectured = Fraction((len(neighbors) + 1) ** 2, 4)
    if k >= 5:
        assert lower > conjectured
    print(
        f"incidence k={k}: order={len(neighbors)}, induced-P5 triples={checked}, "
        f"p(G) in [{half},{half + k}], p(complement)={r}, "
        f"product >= {lower}, conjectured <= {conjectured}, "
        f"strict violation from lower bound={lower > conjectured}"
    )


def main():
    for n in range(1, 7):
        clique = graph(n, combinations(range(n), 2))
        assert by_definition(clique) == n and by_definition(complement(clique)) == 1
    print("source Lemma 3.5: complete and edgeless graphs of orders 1 through 6 verified")

    for k in range(3, 7):
        check_incidence(k)

    for k in range(5, 101):
        r = 2**k
        t = Fraction(2 * k + 1, r)
        next_t = Fraction(2 * k + 3, 2 * r)
        assert next_t <= Fraction(3, 4) * t and t <= Fraction(11, 32)
        assert Fraction((r + 2 * k + 1) ** 2, 4) < Fraction(r * r, 2)
    print("exact arithmetic: ratio recurrence and strict bound checked for k=5 through 100")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
