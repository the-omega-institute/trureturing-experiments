#!/usr/bin/env python3
"""Exact Möbius inversion on multiset partitions of two prime types."""

from functools import cache
from itertools import combinations, product
from math import prod


@cache
def vector_partitions(total, least=(0, 1)):
    """List all sorted multisets of nonzero pairs that sum to total."""
    if total == (0, 0):
        return ((),)
    first_total, second_total = total
    partitions = []
    for a, b in product(range(first_total + 1), range(second_total + 1)):
        pair = (a, b)
        if pair < least:
            continue
        remainder = (first_total - a, second_total - b)
        partitions.extend((pair,) + tail for tail in vector_partitions(remainder, pair))
    return tuple(partitions)


def merge_covers(items, combine):
    """One upper cover for each distinct result of merging two occurrences."""
    covers = set()
    for i, j in combinations(range(len(items)), 2):
        merged = combine(items[i], items[j])
        remaining = items[:i] + items[i + 1:j] + items[j + 1:]
        covers.add(tuple(sorted(remaining + (merged,))))
    return covers


def mobius_values(elements, covers):
    """Compute mu(x, top) and independently mu(bottom, x) on the same order."""
    elements = sorted(elements, key=lambda x: (len(x), x))
    assert len(set(elements)) == len(elements)
    indices = {x: i for i, x in enumerate(elements)}
    assert len(elements[0]) == 1
    assert len(elements) == 1 or len(elements[-2]) < len(elements[-1])
    upper_sets = []
    lower_sets = [set() for _ in elements]
    to_top = []
    for index, element in enumerate(elements):
        upper = set()
        for cover in covers(element):
            cover_index = indices[cover]
            assert cover_index < index
            upper.add(cover_index)
            upper.update(upper_sets[cover_index])
        upper_sets.append(upper)
        for upper_index in upper:
            lower_sets[upper_index].add(index)
        to_top.append(1 if index == 0 else -sum(to_top[i] for i in upper))

    assert len(upper_sets[-1]) == len(elements) - 1
    from_bottom = [0] * len(elements)
    for index in reversed(range(len(elements))):
        from_bottom[index] = (
            1 if index == len(elements) - 1
            else -sum(from_bottom[i] for i in lower_sets[index])
        )
    assert to_top[-1] == from_bottom[0]
    return dict(zip(elements, to_top, strict=True)), dict(zip(elements, from_bottom, strict=True))


def add_pair(x, y):
    return x[0] + y[0], x[1] + y[1]


def check_prime_pair(exponents, vector_to_top, p, q):
    """Check the relabelled integer cover graph and its Möbius values."""
    total = p ** exponents[0] * q ** exponents[1]

    def relabel(factors):
        return tuple(sorted(p ** a * q ** b for a, b in factors))

    integer_elements = {relabel(x) for x in vector_to_top}
    assert len(integer_elements) == len(vector_to_top)
    for factors in vector_to_top:
        integers = relabel(factors)
        assert min(integers) > 1 and prod(integers) == total
        assert {relabel(y) for y in merge_covers(factors, add_pair)} == merge_covers(
            integers, lambda x, y: x * y
        )
    integer_to_top, _ = mobius_values(integer_elements, lambda x: merge_covers(x, lambda a, b: a * b))
    assert all(integer_to_top[relabel(factors)] == value for factors, value in vector_to_top.items())
    bottom = tuple(sorted((p,) * exponents[0] + (q,) * exponents[1]))
    assert integer_to_top[bottom] == 344
    print(f"Relabelling p={p}, q={q}: n={total}, mu(bottom, top)={integer_to_top[bottom]}")


def main():
    for exponents, expected in [((2, 2), -1), ((3, 3), -3), ((4, 4), -8)]:
        to_top, _ = mobius_values(vector_partitions(exponents), lambda x: merge_covers(x, add_pair))
        bottom = ((0, 1),) * exponents[1] + ((1, 0),) * exponents[0]
        assert to_top[bottom] == expected
        print(f"Source table {exponents}: mu(bottom, top)={to_top[bottom]}")

    total = (7, 7)
    to_top, from_bottom = mobius_values(vector_partitions(total), lambda x: merge_covers(x, add_pair))
    bottom = ((0, 1),) * 7 + ((1, 0),) * 7
    top = (total,)
    counts = [sum(len(x) == k for x in from_bottom) for k in range(1, 15)]
    rank_sums = [sum(value for x, value in from_bottom.items() if len(x) == k) for k in range(1, 15)]
    assert counts == [1, 31, 192, 479, 664, 629, 461, 284, 148, 68, 27, 10, 3, 1]
    assert rank_sums == [344, -565, -35, 687, -771, 535, -287, 132, -57, 24, -10, 5, -3, 1]
    assert len(to_top) == 2998 and to_top[bottom] == from_bottom[top] == 344
    assert sum(rank_sums[1:]) == -344
    assert (-1) ** (sum(total) + 1) * to_top[bottom] == -344
    print(f"(7, 7): {len(to_top)} factorizations, mu(bottom, top)={to_top[bottom]}")
    print(f"Counts by number of factors: {counts}")
    print(f"Sum of mu(bottom, x) by number of factors: {rank_sums}")
    for p, q in [(2, 3), (3, 5), (5, 11)]:
        check_prime_pair(total, to_top, p, q)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
