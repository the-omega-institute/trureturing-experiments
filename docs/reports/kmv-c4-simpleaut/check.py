#!/usr/bin/env python3
"""Exact finite checks for KMV, arXiv:2610.22248v1, Conjecture 4.

Only Python's standard library is required. Graphs include loops. A support
permutation is necessary for an invertible structure matrix and sufficient
for some choice of nonzero complex coefficients to be invertible. We count
the entire graph group, an upper bound for the actual permutation image.
These finite checks do not prove the unbounded conjecture.
"""

from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd


def determinant(matrix):
    """Gaussian elimination over Q, independent of support combinatorics."""
    a = [[Fraction(x) for x in row] for row in matrix]
    n, result = len(a), Fraction(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        value = a[j][j]
        result *= value
        for i in range(j + 1, n):
            ratio = a[i][j] / value
            for k in range(j + 1, n):
                a[i][k] -= ratio * a[j][k]
    return result


def rows_of(mask, n):
    return tuple((mask >> (n * i)) & ((1 << n) - 1) for i in range(n))


def distances(rows):
    n = len(rows)
    ds = [-1] * n
    ds[0] = 0
    queue = deque([0])
    while queue:
        i = queue.popleft()
        for j in range(n):
            if (rows[i] >> j) & 1 and ds[j] < 0:
                ds[j] = ds[i] + 1
                queue.append(j)
    return ds


def strongly_connected(rows):
    n = len(rows)
    reach = [rows[i] | (1 << i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if (reach[i] >> k) & 1:
                reach[i] |= reach[k]
    return all(row == (1 << n) - 1 for row in reach)


def cycle_masks(n):
    """All directed simple cycles, with a distinguished smallest vertex."""
    out = []
    for length in range(1, n + 1):
        for vertices in combinations(range(n), length):
            for rest in permutations(vertices[1:]):
                cycle = (vertices[0],) + rest
                mask = sum(1 << (n * cycle[i] + cycle[(i + 1) % length])
                           for i in range(length))
                out.append((mask, length))
    return out


def diagonal_order(rows, ds):
    """Eliminate along actual root-to-vertex paths, then impose every edge.

    Write t_i = z^(2^d_i). The remaining Laurent monomial relations are
    z^(2^(d_i+1) - 2^d_j) = 1. Their integer gcd counts complex solutions.
    No period formula is used by this computation.
    """
    order = 0
    for i, row in enumerate(rows):
        for j in range(len(rows)):
            if (row >> j) & 1:
                order = gcd(order, abs((1 << (ds[i] + 1)) - (1 << ds[j])))
    assert order > 0
    return order


def lattice_order(rows):
    """Independent integer presentation: gcd of all maximal minors.

    The relation row for i -> j is 2 e_i - e_j. The index of this full-rank
    lattice counts its characters in C^*. Used exhaustively through n=3.
    """
    n = len(rows)
    relations = []
    for i, row in enumerate(rows):
        for j in range(n):
            if (row >> j) & 1:
                v = [0] * n
                v[i] += 2
                v[j] -= 1
                relations.append(v)
    index = 0
    for minor in combinations(relations, n):
        value = determinant(minor)
        assert value.denominator == 1
        index = gcd(index, abs(value.numerator))
    assert index > 0
    return index


def graph_data(n):
    perms = list(permutations(range(n)))
    supports = [sum(1 << (n * i + p[i]) for i in range(n)) for p in perms]
    actions = [[sum(((bits >> j) & 1) << p[j] for j in range(n))
                for bits in range(1 << n)] for p in perms]
    return perms, supports, actions


def check_graphs():
    expected = {1: 1, 2: 6, 3: 21, 4: 60}
    for n in range(1, 5):
        perms, supports, actions = graph_data(n)
        cycles = cycle_masks(n)
        count, maximum, lattice_checks = 0, 0, 0
        periods = Counter()
        for mask in range(1 << (n * n)):
            rows = rows_of(mask, n)
            if not strongly_connected(rows):
                continue
            if not any(mask & p == p for p in supports):
                continue
            count += 1
            ds = distances(rows)
            assert min(ds) >= 0
            h = 0
            for c, length in cycles:
                if mask & c == c:
                    h = gcd(h, length)
            assert h >= 1
            classes = [d % h for d in ds]
            sizes = Counter(classes)
            assert set(sizes) == set(range(h))
            assert len(set(sizes.values())) == 1
            k = sizes[0]
            assert n == h * k
            for i in range(n):
                for j in range(n):
                    if (rows[i] >> j) & 1:
                        assert (classes[j] - classes[i] - 1) % h == 0
            diagonal = diagonal_order(rows, ds)
            assert diagonal == (1 << h) - 1
            if n <= 3:
                assert lattice_order(rows) == diagonal
                lattice_checks += 1
            aut = 0
            for p, action in zip(perms, actions):
                if all(action[rows[i]] == rows[p[i]] for i in range(n)):
                    aut += 1
                    rotation = classes[p[0]]
                    assert all((classes[p[i]] - classes[i]) % h == rotation
                               for i in range(n))
            assert aut <= h * factorial(k) ** h
            joint = diagonal * aut
            assert joint <= max(factorial(n), n * ((1 << n) - 1))
            maximum = max(maximum, joint)
            periods[h] += 1
        assert maximum == expected[n]
        print(f"FIDELITY source enumeration: n={n}, admissible_supports={count}, "
              f"max_diagonal_times_graph_aut={maximum}, "
              f"period_counts={dict(sorted(periods.items()))}, "
              f"independent_lattice_checks={lattice_checks}")


def check_attainment():
    orders = []
    for n in range(1, 9):
        cycle = [[int(j == (i + 1) % n) for j in range(n)] for i in range(n)]
        complete = [[1 + int(i == j) for j in range(n)] for i in range(n)]
        assert abs(determinant(cycle)) == 1
        assert determinant(complete) == n + 1
        assert strongly_connected([sum(x << j for j, x in enumerate(row))
                                   for row in cycle])
        assert all(x != 0 for row in complete for x in row)
        cycle_perms = complete_perms = 0
        for p in permutations(range(n)):
            cycle_perms += all(cycle[i][j] == cycle[p[i]][p[j]]
                               for i in range(n) for j in range(n))
            complete_perms += all(complete[i][j] == complete[p[i]][p[j]]
                                  for i in range(n) for j in range(n))
        assert cycle_perms == n
        assert complete_perms == factorial(n)
        m = (1 << n) - 1
        assignments = set()
        for a in range(m):
            phases = tuple((a * (1 << i)) % m for i in range(n))
            assert all((2 * phases[i] - phases[(i + 1) % n]) % m == 0
                       for i in range(n))
            assignments.add(phases)
        assert len(assignments) == m
        # Nonzero loops force t_i = t_i^2, hence t_i = 1 in C^*.
        # The actual lift condition then is coefficient preservation.
        assert cycle_perms * len(assignments) == n * m
        orders.append((n, cycle_perms * m, complete_perms))
    print("ATTAINMENT exact (n, cyclic_full_order, complete_full_order):", orders)


def check_arithmetic():
    count = 0
    for h in range(2, 251):
        for k in range(2, 500 // h + 1):
            lhs = h * ((1 << h) - 1) * factorial(k) ** h
            assert lhs <= factorial(h * k)
            count += 1
    for a in range(1, 65):
        for b in range(1, 65):
            assert gcd((1 << a) - 1, (1 << b) - 1) == (1 << gcd(a, b)) - 1
    print(f"ARITHMETIC factorial comparison: {count} pairs h,k>=2 with hk<=500")
    print("ARITHMETIC Mersenne gcd identity: all 4096 pairs 1<=a,b<=64")


def check_controls():
    # Strong connectivity alone does not imply equal period-class sizes.
    rows = (0b110, 0b001, 0b001)
    assert strongly_connected(rows)
    _, supports, _ = graph_data(3)
    mask = sum(row << (3 * i) for i, row in enumerate(rows))
    assert not any(mask & p == p for p in supports)
    assert Counter(d % 2 for d in distances(rows)) == {0: 1, 1: 2}
    # An actual algebra image can be a proper subgroup of its graph group.
    matrix = [[1, 2], [3, 1]]
    assert determinant(matrix) == -5
    lifts = sum(all(matrix[i][j] == matrix[p[i]][p[j]]
                    for i in range(2) for j in range(2))
                for p in permutations(range(2)))
    assert lifts == 1
    print("CONTROLS singular unequal-period fibers and non-lifting graph symmetry: OK")


if __name__ == "__main__":
    check_graphs()
    check_attainment()
    check_arithmetic()
    check_controls()
    print("ALL_CHECKS_OK")
