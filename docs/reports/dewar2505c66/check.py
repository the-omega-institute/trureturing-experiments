"""Exact checks of the labelled squared-distance map and the bridge bound."""

from fractions import Fraction
from itertools import combinations, product
from math import comb, prod


def squared_distance(left, right):
    return sum((x - y) ** 2 for x, y in zip(left, right, strict=True))


def measurement_image(points, vertex_count, edges):
    distances = {(p, q): squared_distance(p, q) for p, q in product(points, repeat=2)}
    return {
        tuple(distances[placement[u], placement[v]] for u, v in edges)
        for placement in product(points, repeat=vertex_count)
    }


def components(vertex_count, edges):
    neighbors = [set() for _ in range(vertex_count)]
    for u, v in edges:
        neighbors[u].add(v)
        neighbors[v].add(u)
    pending = set(range(vertex_count))
    blocks = set()
    while pending:
        root = pending.pop()
        block = {root}
        frontier = [root]
        while frontier:
            for v in neighbors[frontier.pop()] & pending:
                pending.remove(v)
                block.add(v)
                frontier.append(v)
        blocks.add(frozenset(block))
    return blocks


def triangle_tree(k, tree_edges):
    blocks = [range(3 * i, 3 * i + 3) for i in range(k)]
    internal_edges = tuple(edge for block in blocks for edge in combinations(block, 2))
    joining_edges = tuple((3 * u, 3 * v) for u, v in tree_edges)
    return internal_edges + joining_edges


def check_triangle_fidelity():
    # The source states tight quadratic growth for rigid graphs on equally spaced points.
    edges = ((0, 1), (0, 2), (1, 2))
    for n in range(1, 17):
        image = measurement_image([(i, 0) for i in range(n)], 3, edges)
        # At diameter d > 0, pin the minimum to zero. There are 6d labelled
        # triples attaining both endpoints, paired by reflecting the line.
        assert len(image) == 1 + 3 * comb(n, 2)
    print("fidelity: K3, n=1..16, each image has 1 + 3*binom(n,2) tuples")
    square = tuple(product(range(2), repeat=2))
    assert len(measurement_image(square, 3, edges)) == 10


def check_triangle_tree_structure():
    for k in range(2, 6):
        count = 0
        for tree in combinations(tuple(combinations(range(k), 2)), k - 1):
            if len(components(k, tree)) != 1:
                continue
            count += 1
            edges = triangle_tree(k, tree)
            degrees = [sum(v in edge for edge in edges) for v in range(3 * k)]
            assert min(degrees) == 2
            assert len(components(3 * k, edges)) == 1
            bridges = {
                edge for edge in edges
                if len(components(3 * k, (e for e in edges if e != edge))) == 2
            }
            assert bridges == {(3 * u, 3 * v) for u, v in tree}
            assert components(3 * k, (e for e in edges if e not in bridges)) == {
                frozenset(range(3 * i, 3 * i + 3)) for i in range(k)
            }
        assert count == k ** (k - 2)
        print(f"graph family: k={k}, all {count} labelled trees, {k-1} bridges each")


def check_grid_images():
    triangle_edges = ((0, 1), (0, 2), (1, 2))
    for m, k in [(2, 2), (2, 3), (3, 2)]:
        points = tuple(product(range(m), repeat=2))
        triangle_image = measurement_image(points, 3, triangle_edges)
        grid_distances = {squared_distance(p, q) for p, q in product(points, repeat=2)}
        differences = range(1 - m, m)
        assert grid_distances == {x * x + y * y for x, y in product(differences, repeat=2)}
        edges = triangle_tree(k, tuple((i, i + 1) for i in range(k - 1)))
        image = measurement_image(points, 3 * k, edges)
        for lengths in image:
            assert all(lengths[i : i + 3] in triangle_image for i in range(0, 3 * k, 3))
            assert all(distance in grid_distances for distance in lengths[3 * k :])
        component_bound = len(triangle_image) ** k * len(grid_distances) ** (k - 1)
        normalization_bound = (2 * m - 1) ** (4 * k) * len(grid_distances) ** (k - 1)
        baseline = len(points) ** (3 * k - 1)
        assert len(image) <= component_bound <= normalization_bound
        assert len(image) < baseline
        print(
            f"grid m={m}, k={k}: image={len(image)}, component_bound={component_bound}, "
            f"normalization_bound={normalization_bound}, image/N^(v-1)={Fraction(len(image), baseline)}"
        )


def check_residue_bound():
    primes = (3, 7, 11)
    modulus = prod(p * p for p in primes)
    ratio = prod(Fraction(p * p - p + 1, p * p) for p in primes)
    for p in primes:
        residues = {x * x % (p * p) for x in range(p * p)}
        sums = {(x + y) % (p * p) for x, y in product(residues, repeat=2)}
        assert sums == set(range(p * p)) - set(range(p, p * p, p))
    allowed = {n for n in range(modulus) if all(n % p != 0 or n % (p*p) == 0 for p in primes)}
    assert len(allowed) == ratio * modulus
    for m in (2, 3, 4, 10, 50, 100, 200):
        distances = {x * x + y * y for x, y in product(range(m), repeat=2)}
        assert all(distance % modulus in allowed for distance in distances)
        assert len(distances) <= ratio * (2 * (m - 1) ** 2 + 1) + modulus
    print(f"CRT: primes={primes}, modulus={modulus}, allowed={len(allowed)}, density={ratio}")


if __name__ == "__main__":
    check_triangle_fidelity()
    check_triangle_tree_structure()
    check_grid_images()
    check_residue_bound()
    print("ALL_CHECKS_OK")
