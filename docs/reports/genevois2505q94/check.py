"""Exact finite checks for Genevois, arXiv:2505.09194v1, Question 9.4."""

from collections import deque
from itertools import combinations, product


def subsets(items):
    items = tuple(items)
    for size in range(len(items) + 1):
        yield from map(frozenset, combinations(items, size))


def closure(generators, multiply, identity):
    seen = {identity}
    pending = deque([identity])
    while pending:
        element = pending.popleft()
        for generator in generators:
            neighbour = multiply(element, generator)
            if neighbour not in seen:
                seen.add(neighbour)
                pending.append(neighbour)
    return frozenset(seen)


def conjugate(subgroup, x, inverse, multiply):
    return frozenset(multiply(multiply(x, b), inverse(x)) for b in subgroup)


def check_chain():
    identity = (0, 0)

    def multiply(a, b):
        i, j = a
        k, l = b
        return ((i + pow(5, j, 16) * k) % 16, (j + l) % 4)

    def inverse(a):
        i, j = a
        return ((-pow(13, j, 16) * i) % 16, -j % 4)

    group = frozenset(product(range(16), range(4)))
    lower = frozenset((i, 0) for i in range(16))
    upper = frozenset((0, j) for j in range(4))
    assert pow(5, 4, 16) == pow(13, 4, 16) == 1 and 5 * 13 % 16 == 1
    assert closure(((1, 0), (0, 1)), multiply, identity) == group
    for a, b, c in product(group, repeat=3):
        assert multiply(multiply(a, b), c) == multiply(a, multiply(b, c))
    for a in group:
        assert multiply(a, inverse(a)) == multiply(inverse(a), a) == identity
    for a, b in product(lower, upper):
        twisted = (pow(13, b[1], 16) * a[0] % 16, 0)
        assert multiply(a, b) == multiply(b, twisted)
    print("Example 2.6: two-factor semidirect model, 64 elements, 64 twist relations")

    factors = (lower, upper)
    standards = {
        indexes: closure(set().union(*(factors[i] for i in indexes)), multiply, identity)
        for indexes in subsets(range(len(factors)))
    }
    for r, s in product(standards, repeat=2):
        assert standards[r] & standards[s] == standards[r & s]
    print("Corollary 3.7: standard intersection formula, all 16 stable pairs")

    x = (2, 0)
    shifted = conjugate(upper, x, inverse, multiply)
    expected = frozenset({(0, 0), (0, 2)})
    assert shifted == {(0, 0), (8, 1), (0, 2), (8, 3)}
    assert upper & shifted == expected
    parabolics = {
        conjugate(standard, g, inverse, multiply)
        for standard, g in product(standards.values(), group)
    }
    assert {len(p) for p in parabolics} == {1, 4, 16, 64}
    assert expected not in parabolics
    print(f"Question 9.4: intersection {sorted(expected)}, parabolic sizes [1, 4, 16, 64]")

    # In this Cayley graph an edge is multiplication on the right by a factor element.
    generators = (lower | upper) - {identity}
    coset = {multiply(x, b) for b in upper}
    projected = set()
    for v in coset:
        neighbours = {multiply(v, generator) for generator in generators}
        closest = neighbours & upper
        assert len(closest) == 1
        gate = next(iter(closest))
        assert {multiply(gate, generator) for generator in upper - {identity}} == upper - {gate}
        assert not (v in upper or neighbours & (upper - {gate}))
        projected.add(gate)
    assert projected == upper
    projection_stabilizer = {g for g in group if {multiply(g, p) for p in projected} == projected}
    both_stabilizer = projection_stabilizer & {
        g for g in group if {multiply(g, v) for v in coset} == coset
    }
    assert projection_stabilizer == upper and both_stabilizer == expected
    print("Question 9.4 geometric sketch: stabilizer sizes 4 and 2")


def check_lamps(n):
    positions = range(n)
    mask = (1 << n) - 1
    identity = (0, 0)

    def translate(lamps, shift):
        shift %= n
        return ((lamps << shift) | (lamps >> (n - shift))) & mask

    def multiply(a, b):
        lamps, shift = a
        other, step = b
        return (lamps ^ translate(other, shift), (shift + step) % n)

    def inverse(a):
        lamps, shift = a
        return (translate(lamps, -shift), -shift % n)

    upper = {(0, j) for j in positions}
    lamps = [(1 << i, 0) for i in positions]
    for lamp in lamps:
        assert multiply(lamp, lamp) == identity
    for a, b in product(lamps, repeat=2):
        assert multiply(a, b) == multiply(b, a)
    for i, j in product(positions, repeat=2):
        assert multiply(lamps[i], (0, j)) == multiply((0, j), lamps[(i - j) % n])
    assert len({tuple(translate(a, j) for a, _ in lamps) for j in positions}) == n
    assert len(closure([*lamps, (0, 1)], multiply, identity)) == (1 << n) * n

    orders = []
    for k in range(1, n + 1):
        if n % k:
            continue
        subgroup = set(range(0, n, n // k))
        indicator = sum(1 << i for i in subgroup)
        shifted = conjugate(upper, (indicator, 0), inverse, multiply)
        assert upper & shifted == {(0, j) for j in subgroup}
        orders.append(k)

    # The top-coordinate image of a generated subgroup is generated by the
    # top coordinates of its factor generators. Conjugation in C_n is trivial.
    projection_sizes = set()
    for indexes in subsets(range(n + 1)):
        projected_generators = {int(i == n) for i in indexes}
        image = closure(projected_generators, lambda a, b: (a + b) % n, 0)
        projection_sizes.add(len(image))
    assert projection_sizes == {1, n}
    assert all(k not in projection_sizes for k in orders if 1 < k < n)
    print(f"Lamp family C_{n}: realized orders {orders}, parabolic image sizes [1, {n}]")


if __name__ == "__main__":
    check_chain()
    for size in range(2, 9):
        check_lamps(size)
    print("ALL_CHECKS_OK")
