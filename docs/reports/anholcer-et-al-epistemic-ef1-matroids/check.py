#!/usr/bin/env python3
"""Exact checks for arXiv:2606.11494v1, Conjecture 4.

Only Python integers and fractions are used. The finite checks are evidence for
the preregistration; they do not constitute a Lean proof of the unbounded family.
"""

from fractions import Fraction as Q
from itertools import combinations, product


def value(bundle, weights):
    return sum((weights[i] for i in bundle), Q(0))


def goods_ef1(allocation, valuations):
    return all(
        not other
        or any(value(other - {g}, w) <= value(own, w) for g in other)
        for own, w in zip(allocation, valuations)
        for other in allocation
    )


def chores_ef1(allocation, valuations):
    return all(
        i == j or not own
        or any(value(own - {c}, w) <= value(other, w) for c in own)
        for i, (own, w) in enumerate(zip(allocation, valuations))
        for j, other in enumerate(allocation)
    )


def forest(edges):
    parent = {}

    def root(v):
        parent.setdefault(v, v)
        while parent[v] != v:
            v = parent[v]
        return v

    for u, v in edges:
        a, b = root(u), root(v)
        if a == b:
            return False
        parent[a] = b
    return True


def fidelity():
    # Section 4, "EXAMPLE (illustration of both theorems)".
    edges = [(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (4, 6)]
    bundles = [frozenset({0, 1, 2}), frozenset({3, 4}), frozenset({5, 6})]
    weights = [(3, 3, 3, 1, 1, 1, 1),
               (1, 1, 1, 4, 4, 1, 1),
               (1, 1, 1, 1, 1, 1, 1)]
    assert set.union(*(set(s) for s in bundles)) == set(range(7))
    assert all(not (a & b) for a, b in combinations(bundles, 2))
    assert all(forest([edges[e] for e in b]) for b in bundles)
    matrix = [[value(b, w) for b in bundles] for w in weights]
    assert matrix == [[9, 2, 2], [3, 8, 2], [3, 2, 2]]
    assert goods_ef1(bundles, weights)
    # Every cyclic length-five block is a spanning tree, as the source states.
    assert all(forest([edges[(start + k) % 7] for k in range(5)])
               for start in range(7))
    left, right = frozenset({0, 1, 2}), frozenset({3, 4, 5, 6})

    def monotone_value(s):
        return sum(3 in edges[e] for e in s) ** 2 + sum(4 in edges[e] for e in s) ** 2

    assert monotone_value(left) == monotone_value(right) == 2
    print("FIDELITY: source seven-edge example values [[9,2,2],[3,8,2],[3,2,2]]; goods-EF1")
    print("FIDELITY: seven cyclic five-edge blocks are trees; source nonlinear values 2,2")


def mask_set(mask, size):
    return frozenset(i for i in range(size) if mask & (1 << i))


def independent(mask, blocks):
    return all((mask & block).bit_count() <= 1 for block in blocks)


def matroid_axioms(indep, size):
    # Exhaustive finite independence axioms, including every exchange pair.
    family = set(indep)
    assert 0 in family
    for a in indep:
        for e in range(size):
            assert a & ~(1 << e) in family
    for a in indep:
        for b in indep:
            if a.bit_count() > b.bit_count():
                assert any((a & ~b) & (1 << e) and b | (1 << e) in family
                           for e in range(size))


def cycle(r):
    size = 2 * r
    full = (1 << size) - 1
    blocks1 = [(1 << i) | (1 << (i + 1)) for i in range(0, size, 2)]
    blocks2 = [(1 << i) | (1 << ((i + 1) % size)) for i in range(1, size, 2)]
    families = [[s for s in range(full + 1) if independent(s, blocks)]
                for blocks in (blocks1, blocks2)]
    odd = sum(1 << i for i in range(0, size, 2))  # one-based odd objects
    even = full ^ odd
    for blocks, family in zip((blocks1, blocks2), families):
        assert max(s.bit_count() for s in family) == r
        assert all(independent(1 << e, blocks) for e in range(size))
        assert not independent(full, blocks)
        assert independent(odd, blocks) and independent(even, blocks)
        # The preceding two checks give chromatic number exactly two.
        if r <= 4:
            matroid_axioms(family, size)
    # Enumerate every complete labeled allocation, without assuming saturation.
    feasible = [s for s in range(full + 1)
                if independent(s, blocks1) and independent(full ^ s, blocks2)]
    assert set(feasible) == {odd, even}
    allocations = [(mask_set(s, size), mask_set(full ^ s, size)) for s in feasible]
    for b, epsilon in product((Q(1), Q(2, 3), Q(5, 2)), (Q(1, 10), Q(1), Q(7, 3))):
        a = r * b / (r - 1) + epsilon
        weights = tuple(a if i % 2 == 0 else b for i in range(size))
        assert (r - 1) * a > r * b
        for allocation in allocations:
            assert not goods_ef1(allocation, [weights, weights])
            assert not chores_ef1(allocation, [weights, weights])
    # Equality and the lower side of the threshold are controls, not refutations.
    for a in (Q(1), Q(r, r - 1)):
        weights = tuple(a if i % 2 == 0 else Q(1) for i in range(size))
        assert all(goods_ef1(x, [weights, weights]) for x in allocations)
        assert all(chores_ef1(x, [weights, weights]) for x in allocations)
    print(f"CYCLE r={r}: rank={r}, loopless, chi=2 each; exactly 2 feasible allocations; "
          "9 rational parameter pairs fail goods/chores EF1; threshold controls pass")
    return allocations


def four_items():
    allocations = cycle(2)
    weights = (3, 1, 3, 1)
    for allocation in allocations:
        assert sorted(value(s, weights) for s in allocation) == [2, 6]
        high = next(s for s in allocation if value(s, weights) == 6)
        assert all(value(high - {e}, weights) == 3 for e in high)
    print("WITNESS E={1,2,3,4}: values 6 and 2; every high-bundle deletion leaves 3 > 2")
    # Every corner of an independent eight-coordinate perturbation box.
    # Each of the two agents has its own perturbed singleton values.
    for signs in product((-1, 1), repeat=8):
        valuations = [tuple(Q(weights[e]) + Q(signs[4 * i + e], 10)
                            for e in range(4)) for i in range(2)]
        assert all(x > 0 for w in valuations for x in w)
        for x in allocations:
            assert not goods_ef1(x, valuations)
            assert not chores_ef1(x, valuations)
    print("ROBUSTNESS: all 256 corners of the positive 8-coordinate radius-1/10 box fail both EF1 notions")


def main():
    fidelity()
    four_items()
    for r in range(3, 9):
        cycle(r)
    print("MATROID AXIOMS: exhaustive heredity and exchange for each matroid at r=2,3,4")
    print("SCOPE: finite exact checks for r=2..8; infinite-family argument is in the preregistration")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
