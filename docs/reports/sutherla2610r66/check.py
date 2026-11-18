#!/usr/bin/env python3
"""Exact checks for Sutherland, arXiv:2610.01049v4, Remark 6.6."""

from collections import deque
from itertools import product

# Signed indices encode 1, i, j, k. The Boolean factor uses bit vectors.
HAMILTON = (
    (1, 2, 3, 4),
    (2, -1, 4, -3),
    (3, -4, -1, 2),
    (4, 3, -2, -1),
)
QUATERNIONS = (1, -1, 2, -2, 3, -3, 4, -4)
IDENTITY = (1, 0)
CENTRAL_SQUARE = (-1, 0)


def multiply(x, y):
    q, u = x
    r, v = y
    sign = -1 if (q < 0) != (r < 0) else 1
    return sign * HAMILTON[abs(q) - 1][abs(r) - 1], u ^ v


def inverse(x):
    q, u = x
    return (q if abs(q) == 1 else -q), u


def right_difference(x, y):
    return multiply(x, inverse(y))


def inverse_pair(x):
    return frozenset((x, inverse(x)))


def generated_subgroup(generators):
    seen = {IDENTITY}
    frontier = deque(seen)
    while frontier:
        x = frontier.popleft()
        for s in generators:
            y = multiply(x, s)
            if y not in seen:
                seen.add(y)
                frontier.append(y)
    return seen


def check_codegree(vertices, connection):
    """Compare graph codegrees, translated intersections, and ordered differences."""
    neighbors = {
        x: {y for y in vertices if right_difference(x, y) in connection}
        for x in vertices
    }
    for x, y in product(vertices, repeat=2):
        if x == y:
            continue
        g = right_difference(y, x)
        common_neighbors = len(neighbors[x] & neighbors[y])
        translated_intersection = len(connection & {multiply(s, g) for s in connection})
        representations = sum(
            multiply(inverse(u), v) == g for u, v in product(connection, repeat=2)
        )
        assert common_neighbors == translated_intersection == representations
    print(f"fidelity codegree pairs={len(vertices) * (len(vertices) - 1)} (Section 1.1)")


def rooted_xi_size(vertices, labels):
    """Count rooted permutations satisfying the complete difference coloring."""
    assigned = {IDENTITY: IDENTITY}
    used = {IDENTITY}

    def count(index):
        if index == len(vertices):
            return 1
        x = vertices[index]
        total = 0
        for image in inverse_pair(x):
            if image in used:
                continue
            if any(labels[x, y] != labels[image, target] for y, target in assigned.items()):
                continue
            assigned[x] = image
            used.add(image)
            total += count(index + 1)
            used.remove(image)
            del assigned[x]
        return total

    return count(1)


def check(rank):
    vertices = list(product(QUATERNIONS, range(1 << rank)))
    labels = {(x, y): inverse_pair(right_difference(x, y)) for x, y in product(vertices, repeat=2)}
    assert rooted_xi_size(vertices, labels) == 8
    print(f"fidelity rank={rank} vertices={len(vertices)} [Xi:R]=8 (Theorem 1.1)")

    if rank == 0:
        return

    steps = {s for s in vertices if multiply(s, s) != IDENTITY}
    assert all(multiply(s, s) == CENTRAL_SQUARE for s in steps)
    assert {inverse(s) for s in steps} == steps
    assert {multiply(CENTRAL_SQUARE, s) for s in steps} == steps
    assert multiply((2, 0), (3, 0)) != multiply((3, 0), (2, 0))
    assert (2, 0) in steps and (3, 0) in steps

    generators = {(q, 0) for q in (2, -2, 3, -3)}
    generators |= {(q, 1 << n) for q in (2, -2) for n in range(rank)}
    assert generators <= steps
    assert {inverse(s) for s in generators} == generators
    assert generated_subgroup(generators) == set(vertices)
    check_codegree(vertices, generators)

    swap = {x: x for x in vertices}
    swap[IDENTITY], swap[CENTRAL_SQUARE] = CENTRAL_SQUARE, IDENTITY
    b = (1, 1)
    assert labels[IDENTITY, b] == frozenset((b,))
    assert labels[swap[IDENTITY], swap[b]] == frozenset((multiply(CENTRAL_SQUARE, b),))
    assert labels[IDENTITY, b] != labels[swap[IDENTITY], swap[b]]

    # Each admissible connection set chooses an arbitrary subset of these
    # inverse pairs. Giving every other difference label the same value
    # therefore checks all such connection sets at once.
    step_labels = {inverse_pair(s) for s in steps}
    eligible_labels = {
        edge: label if label in step_labels else None for edge, label in labels.items()
    }
    for x, y in product(vertices, repeat=2):
        assert eligible_labels[x, y] == eligible_labels[swap[x], swap[y]]
    print(
        f"obstruction rank={rank} inverse_pairs={len(step_labels)} "
        f"connection_sets={1 << len(step_labels)} generators={len(generators)}"
    )


def main():
    for a, b, c in product(QUATERNIONS, repeat=3):
        x, y, z = (a, 0), (b, 0), (c, 0)
        assert multiply(multiply(x, y), z) == multiply(x, multiply(y, z))
        assert multiply(x, inverse(x)) == IDENTITY
        assert multiply(IDENTITY, x) == x == multiply(x, IDENTITY)
    for rank in range(5):
        check(rank)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
