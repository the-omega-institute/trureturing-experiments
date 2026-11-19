#!/usr/bin/env python3
"""Check a source lemma and sample the coupled Heuberger matrix family."""

from itertools import combinations, product
from math import gcd


def chromatic_number(adjacency):
    """Return the least number of colors, or None if a loop prevents coloring."""
    if any(vertex in neighbors for vertex, neighbors in enumerate(adjacency)):
        return None

    colors = [-1] * len(adjacency)

    def extend(vertex, count):
        if vertex == len(adjacency):
            return True
        forbidden = {colors[neighbor] for neighbor in adjacency[vertex]}
        for color in range(count):
            if color not in forbidden:
                colors[vertex] = color
                if extend(vertex + 1, count):
                    return True
        colors[vertex] = -1
        return False

    for count in range(1, len(adjacency) + 1):
        if extend(0, count):
            return count
    raise AssertionError("each vertex can be assigned its own color")


def one_row_cayley(first, second):
    """Quotient by both columns via connected orbits in a finite cyclic group."""
    modulus = abs(first or second)
    representatives = {}
    for start in range(modulus):
        if start in representatives:
            continue
        orbit = [start]
        representatives[start] = start
        for vertex in orbit:
            for step in (first, second, -first, -second):
                neighbor = (vertex + step) % modulus
                if neighbor not in representatives:
                    representatives[neighbor] = start
                    orbit.append(neighbor)

    vertices = sorted(set(representatives.values()))
    indices = {v: i for i, v in enumerate(vertices)}
    return [
        {indices[representatives[(v + step) % modulus]] for step in (-1, 1)}
        for v in vertices
    ]


def check_source():
    """Reproduce the one-row gcd lemma in arXiv:2511.03028v1, section 2.1."""
    counts = {None: 0, 2: 0, 3: 0}
    for first, second in product(range(-12, 13), repeat=2):
        if first == second == 0:
            continue
        order = gcd(first, second)
        expected = None if order == 1 else 2 if order % 2 == 0 else 3
        adjacency = one_row_cayley(first, second)
        actual = chromatic_number(adjacency)
        assert actual == expected, (first, second, expected, actual)
        counts[actual] += 1
    print(
        f"Source §2.1: {sum(counts.values())} one-row matrices, "
        f"{counts[None]} looped, {counts[2]} with χ=2, {counts[3]} with χ=3"
    )


def identity(size):
    return [tuple(int(i == j) for j in range(size)) for i in range(size)]


def check_family(rank, blocks):
    tail_size = 4 * blocks
    dimension = rank + tail_size
    basis = identity(dimension)

    columns = [
        tuple(value - multiple * first for value, first in zip(e, basis[0]))
        for multiple, e in enumerate(basis[1:rank], start=2)
    ]
    columns.append((rank + 1,) + (0,) * (rank - 1) + tuple(range(1, tail_size + 1)))
    rows = list(zip(*columns))
    pivots = list(range(1, rank + 1))
    free = [0, *range(rank + 1, dimension)]

    assert [rows[i] for i in pivots] == identity(rank), (rank, blocks)
    assert all(any(row) for row in rows) and len(set(rows)) == dimension
    assert all(gcd(*column) == 1 for column in columns)

    # Eliminate the pivot coordinates with the input columns. The identity
    # minor makes this integral and gives a section using the free coordinates.
    def quotient(vector):
        return tuple(
            vector[i] - sum(vector[pivot] * entry for pivot, entry in zip(pivots, rows[i]))
            for i in free
        )

    assert all(not any(quotient(column)) for column in columns)
    assert [quotient(basis[i]) for i in free] == identity(tail_size)
    generator_images = [quotient(e) for e in basis]
    signed = set(generator_images)
    signed.update(tuple(-entry for entry in image) for image in generator_images)
    assert (0,) * tail_size not in signed and len(signed) == 2 * dimension

    # Check the clique by exact adjacency in the integral quotient.
    vertices = [quotient((t,) + (0,) * (dimension - 1)) for t in range(rank + 1)]
    assert len(set(vertices)) == rank + 1
    for u, v in combinations(vertices, 2):
        assert tuple(x - y for x, y in zip(u, v)) in signed

    # Check that the proposed coloring is a function on cosets. Test every
    # generator and its negative using the coloring on the free coordinates.
    colors = tuple(range(1, rank + 1)) + (1, -1, -1, 1) * blocks
    modulus = rank + 1
    assert all(sum(a * b for a, b in zip(colors, column)) % modulus == 0 for column in columns)
    free_colors = tuple(colors[i] for i in free)
    assert all(sum(a * b for a, b in zip(free_colors, g)) % modulus != 0 for g in signed)
    assert len({sum(a * b for a, b in zip(free_colors, v)) % modulus for v in vertices}) == modulus


def main():
    check_source()
    for rank, blocks in product(range(3, 13), range(1, 9)):
        check_family(rank, blocks)
    print("Family: 80 pairs (r,k), 3 ≤ r ≤ 12 and 1 ≤ k ≤ 8, verified χ=r+1")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
