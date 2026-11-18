"""Exact checks for Slavković–Zhu–Petrović, arXiv:1401.1397v1, Conjecture 3.3."""

from fractions import Fraction
from itertools import combinations
from math import gcd, lcm


def denominator_row(slices):
    """Each slice is a two-row A-by-C table at one B level."""
    row = []
    for counts in slices:
        marginals = tuple(map(sum, counts))
        total = sum(marginals)
        row.append(lcm(*(Fraction(count, total).denominator for count in marginals)))
    return tuple(row)


def realize(weights):
    return tuple(((1, 2), (weight - 2, 2 * weight - 1)) for weight in weights)


def check_realization(weights):
    slices = realize(weights)
    for (a, b), (c, d) in slices:
        assert min(a, b, c, d) > 0
        assert a * d - b * c == 3
    assert denominator_row(slices) == weights


def support_components(weights, reachable, degree):
    """Components of the support intersection graph on a nonzero fiber."""
    unseen = {
        i for i, weight in enumerate(weights)
        if degree >= weight and reachable[degree - weight]
    }
    num_components = 0
    while unseen:
        num_components += 1
        frontier = [unseen.pop()]
        while frontier:
            i = frontier.pop()
            adjacent = {
                j for j in unseen
                if degree >= weights[i] + weights[j]
                and reachable[degree - weights[i] - weights[j]]
            }
            unseen.difference_update(adjacent)
            frontier.extend(adjacent)
    return num_components


def minimum_markov_size(weights):
    """Sum the number of inter-component edges needed in each fiber."""
    assert len(weights) >= 2 and min(weights) > 0
    scale = gcd(*weights)
    weights = tuple(weight // scale for weight in weights)
    reachable = [True]
    consecutive = 1
    smallest = min(weights)
    while consecutive < smallest:
        degree = len(reachable)
        present = any(degree >= w and reachable[degree - w] for w in weights)
        reachable.append(present)
        consecutive = consecutive + 1 if present else 0
    conductor = len(reachable) - smallest
    # Above this cutoff, any two active indices appear together in a factorization.
    cutoff = conductor + max(sum(pair) for pair in combinations(weights, 2))
    reachable.extend([True] * (cutoff - len(reachable)))
    needs = []
    for degree in range(1, cutoff):
        components = support_components(weights, reachable, degree)
        if components > 1:
            needs.append((degree * scale, components - 1))
    return sum(count for _, count in needs), cutoff * scale, needs


def main():
    # The original counts, grouped by B (gender), from arXiv §4.1, Table 1.
    source = (((8, 7), (4, 6)), ((2, 3), (9, 11)))
    assert sum(sum(row) for group in source for row in group) == 50
    weights = denominator_row(source)
    assert weights == (5, 5)
    positive_solutions = [
        (x, y)
        for x in range(1, 11)
        for y in range(1, 11)
        if weights[0] * x + weights[1] * y == 50
    ]
    assert positive_solutions == [(x, 10 - x) for x in range(1, 10)]
    assert minimum_markov_size(weights)[0] == 1
    print("Source §4.1: weights=(5, 5), total=50, positive_solutions=9, Markov size=1")

    for n in range(3, 13):
        weights = tuple(range(n, 2 * n))
        check_realization(weights)
        size, cutoff, needs = minimum_markov_size(weights)
        assert size == n * (n - 1) // 2 > n - 1
        if n == 3:
            assert needs == [(8, 1), (9, 1), (10, 1)]
        print(f"Consecutive n={n}: minimum={size}, cutoff={cutoff}, split_fibers={len(needs)}")

    for h in range(2, 13):
        a, b, c = 2 * h - 1, 2 * h, 2 * h + 1
        weights = (a * b, a * c, b * c, b * c + a)
        check_realization(weights)
        size, cutoff, needs = minimum_markov_size(weights)
        assert size == 4 * h > 3
        print(f"Bresinsky h={h}: minimum={size}, cutoff={cutoff}, split_fibers={len(needs)}")

    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
