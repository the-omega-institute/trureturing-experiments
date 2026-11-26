"""Exact witnesses for the two matrix condensates and the quartic obstruction."""
from collections import Counter
from itertools import combinations, permutations, product

PAIRS = tuple(combinations(range(5), 2))
TRIANGLE = frozenset(combinations(range(3), 2))
STAR = frozenset((0, leaf) for leaf in (1, 2, 3))


def phase_ranks():
    source = {p: sum(i < 3 for i in p) % 3 for p in PAIRS}
    target = {p: int(sum(i < 3 for i in p) == 1) for p in PAIRS}
    s, t = Counter(source.values()), Counter(target.values())
    ranks = [sum(min(n, t[(phase + q) % 3]) for phase, n in s.items())
             for q in range(3)]
    assert ranks == [7, 4, 7]
    assert all((target[p] - source[p]) % 3 == int(p in TRIANGLE) for p in PAIRS)
    assert sum(target.values()) % 3 == 0
    return ranks


def operator_entry(support, upper, lower):
    if upper[0] == upper[1] or tuple(sorted(upper)) not in support:
        return 0
    return int(upper == lower) - int(upper == lower[::-1])


def contractions(left, right):
    # The 24 possible pairings of four upper with four lower fundamental indices.
    values = []
    for perm in permutations(range(4)):
        total = 0
        for indices in product(range(5), repeat=4):
            lower = tuple(indices[i] for i in perm)
            total += (operator_entry(left, indices[:2], lower[:2])
                      * operator_entry(right, indices[2:], lower[2:]))
        values.append(total)
    return values


def verify():
    ranks = phase_ranks()
    all_pairs = frozenset(PAIRS)
    for use_x_left, use_x_right in product((False, True), repeat=2):
        values = []
        for shape in (TRIANGLE, STAR):
            left = shape if use_x_left else all_pairs - shape
            right = shape if use_x_right else all_pairs - shape
            values.append(contractions(left, right))
        assert values[0] == values[1]
    spectra = [sorted((sum(i in edge for edge in shape) for i in range(5)),
                      reverse=True) for shape in (TRIANGLE, STAR)]
    assert spectra == [[2, 2, 2, 0, 0], [3, 1, 1, 1, 0]]
    assert [sum(t * (t - 2)**2 for t in spec) for spec in spectra] == [0, 6]
    assert [sum(t*t for t in spec) for spec in spectra] == [12, 12]
    # Exhaust the coordinate three-planes as an additional check of the sextic.
    zeros = []
    for support in combinations(PAIRS, 3):
        d = [sum(i in edge for edge in support) for i in range(5)]
        if sum(t * (t - 2)**2 for t in d) == 0:
            zeros.append(frozenset(support))
    assert set(zeros) == {frozenset(combinations(c, 2)) for c in combinations(range(5), 3)}
    print("Single-field maximal ranks", ranks)
    print("All 96 quartic coordinate pairings agree. Sextic energies at v=1: 0 and 6.")


if __name__ == "__main__":
    verify()
