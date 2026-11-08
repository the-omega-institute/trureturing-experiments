#!/usr/bin/env python3
"""Exact finite checks for Mehdi-Nezhad, arXiv:2611.08022v1, Problem 7.2.

Python standard library only. Families are enumerated as maximal cliques in
an intersection graph, independently of the proposed threshold construction.
All numerical evidence uses integers or fractions.Fraction.
"""
from fractions import Fraction as Q
from itertools import product
from random import Random
from math import prod

SEED = 261108022


def maximal_intersecting(k):
    vertices = list(range(1, (1 << k) - 1))
    adjacency = [sum(1 << j for j, b in enumerate(vertices)
                     if a != b and a & b) for a in vertices]
    answers = []

    def visit(chosen, possible, excluded):
        if not (possible | excluded):
            answers.append(tuple(vertices[i] for i in range(len(vertices))
                                 if chosen >> i & 1))
            return
        union = possible | excluded
        indices = [i for i in range(len(vertices)) if union >> i & 1]
        pivot = max(indices, key=lambda i: (possible & adjacency[i]).bit_count())
        candidates = possible & ~adjacency[pivot]
        while candidates:
            bit = candidates & -candidates
            i = bit.bit_length() - 1
            visit(chosen | bit, possible & adjacency[i], excluded & adjacency[i])
            possible &= ~bit
            excluded |= bit
            candidates &= ~bit

    visit(0, (1 << len(vertices)) - 1, 0)
    assert len(set(answers)) == len(answers)
    for family in answers:
        assert all(a & b for a in family for b in family)
        assert all(a in family or any(not (a & b) for b in family) for a in vertices)
        assert len(family) == (1 << (k - 1)) - 1
    return answers


def subsets(mask):
    a = mask
    while True:
        yield a
        if a == 0:
            return
        a = (a - 1) & mask


def tprod(t, mask):
    return prod((x for i, x in enumerate(t) if mask >> i & 1), start=Q(1))


def envelope(c, t):
    full = (1 << len(t)) - 1
    low = sum(1 << i for i, x in enumerate(t) if x < 1)
    high = full ^ low
    return (c / 2 * prod((1 + x for i, x in enumerate(t) if high >> i & 1), start=Q(1))
            * sum((min(tprod(t, b), tprod(t, low ^ b)) for b in subsets(low)), Q(0))
            - c * tprod(t, full))


def threshold(t, ell):
    full = (1 << len(t)) - 1
    low = sum(1 << i for i, x in enumerate(t) if x < 1)
    assert low >> ell & 1
    return tuple(a for a in range(1, full)
                 if tprod(t, a & low) < tprod(t, low ^ (a & low))
                 or (tprod(t, a & low) == tprod(t, low ^ (a & low))
                     and a >> ell & 1))


def source_low(c, t):
    # Theorem 4.4: c_source = c * product(t), a_j = 1/t_j.
    full = (1 << len(t)) - 1
    a = tuple(1 / x for x in t)
    return c * tprod(t, full) / 2 * sum(
        (min(tprod(a, b), tprod(a, full ^ b)) for b in range(1, full)), Q(0))


def source_high(c, t):
    i = min(range(len(t)), key=lambda j: t[j])
    return c * (t[i] * prod((1 + x for j, x in enumerate(t) if j != i), start=Q(1))
                - tprod(t, (1 << len(t)) - 1))


def check_vector(families, c, t):
    full = (1 << len(t)) - 1
    weights = [c * tprod(t, a) for a in range(full + 1)]
    costs = [sum((weights[a] for a in f), Q(0)) for f in families]
    true_min = min(costs)
    low_indices = [i for i, x in enumerate(t) if x < 1]
    if low_indices:
        assert true_min == envelope(c, t), (t, true_min, envelope(c, t))
        for ell in low_indices:
            f = threshold(t, ell)
            assert f in families
            assert sum((weights[a] for a in f), Q(0)) == true_min
        # Stronger optional diagnostic: strictly expensive coordinates are ignored
        # on the full cube by every optimizer, not only by our construction.
        for f, cost in zip(families, costs):
            if cost != true_min:
                continue
            cube_family = set(f) | {full}
            for i, x in enumerate(t):
                if x > 1:
                    assert all((a in cube_family) == ((a ^ (1 << i)) in cube_family)
                               for a in range(full + 1))
    if all(x <= 1 for x in t):
        assert true_min == source_low(c, t)
        if low_indices:
            assert envelope(c, t) == source_low(c, t)
    if all(x >= 1 for x in t):
        assert true_min == source_high(c, t)
        # L is empty here: the envelope is valid only on the neutral boundary.
        if min(t) == 1:
            assert envelope(c, t) == source_high(c, t)
    return true_min


def main():
    counts = {2: 2, 3: 4, 4: 12, 5: 81, 6: 2646}
    all_families = {}
    for k, expected in counts.items():
        all_families[k] = maximal_intersecting(k)
        assert len(all_families[k]) == expected
        print(f'ENUMERATION k={k} maximal_intersecting={expected}')
    f3 = all_families[3]
    t = (Q(2),) * 3
    costs = sorted(sum((tprod(t, a) for a in f), Q(0)) for f in f3)
    assert costs == [10, 10, 10, 12]
    pair_bound = Q(1, 2) * sum(min(tprod(t, a), tprod(t, 7 ^ a)) for a in range(1, 7))
    assert pair_bound == 6
    print('FIDELITY Remark 4.5: weights 2^|A|; family costs=10,10,10,12; pair bound=6')
    for k in counts:
        assert check_vector(all_families[k], Q(1), (Q(1),) * k) == (1 << (k - 1)) - 1
    print('FIDELITY Boolean weights: i=2^(k-1)-1 for k=2..6')
    t = (Q(1, 2), Q(1, 2), Q(1, 2), Q(2))
    assert check_vector(all_families[4], Q(8), t) == 19
    f = threshold(t, 0)
    assert set(f) == {3, 5, 6, 7, 11, 13, 14}
    stars = [sum((8 * tprod(t, a) for a in range(1, 15) if a >> i & 1), Q(0))
             for i in range(4)]
    assert stars == [25, 25, 25, 52]
    print('MIXED INSTANCE c=8,t=(1/2,1/2,1/2,2): minimum=19; stars=25,25,25,52')
    grid_count = 0
    for k in [2, 3, 4]:
        for t in product((Q(1, 2), Q(1), Q(2)), repeat=k):
            check_vector(all_families[k], Q(3, 2), t)
            grid_count += 1
    print(f'EXHAUSTIVE PARAMETER GRID k=2..4 over (1/2,1,2): vectors={grid_count}')
    rng = Random(SEED)
    pool = (Q(1, 4), Q(1, 2), Q(2, 3), Q(1), Q(3, 2), Q(2), Q(3), Q(4))
    for k in counts:
        for _ in range(33):
            t = [pool[rng.randrange(3)], pool[rng.randrange(4, len(pool))]]
            t += [rng.choice(pool) for _ in range(k - 2)]
            rng.shuffle(t)
            check_vector(all_families[k], Q(rng.randint(1, 9), rng.randint(1, 5)), tuple(t))
        print(f'MIXED EXACT seed={SEED}, k={k}: vectors=33, mismatches=0')
    # Additional unmixed checks verify both stated source regimes through k=6.
    for k in counts:
        for _ in range(8):
            low = tuple(rng.choice(pool[:4]) for _ in range(k))
            high = tuple(rng.choice(pool[3:]) for _ in range(k))
            check_vector(all_families[k], Q(2), low)
            check_vector(all_families[k], Q(2), high)
    print('UNMIXED EXACT k=2..6: vectors=80; Theorem 4.4 and Proposition 6.5 agree')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
