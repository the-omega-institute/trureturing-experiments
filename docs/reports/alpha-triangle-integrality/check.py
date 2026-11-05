#!/usr/bin/env python3
"""Exact small-case checks for arXiv:1203.4419v1, Conjecture 3.1.

Enumerate labelled Ferrers diagrams by adding every admissible cell, not by
assuming the conjectured polynomial or stabilizer divisibility. No dependencies.
"""
from fractions import Fraction
from itertools import permutations
from math import comb, factorial

MAX_M = 4
SOURCE_C = [
    [1],
    [0, 1, 1],
    [0, 1, 3, 6, 3],
    [0, 1, 7, 20, 46, 45, 15],
    [0, 1, 11, 61, 198, 480, 645, 420, 105],
]


def mu(x):
    return frozenset([(0,) * x] + [tuple(int(i == j) for i in range(x)) for j in range(x)])


def grow(diagrams, x):
    result = set()
    for diagram in diagrams:
        candidates = set()
        for v in diagram:
            for i in range(x):
                candidates.add(tuple(a + int(j == i) for j, a in enumerate(v)))
        for v in candidates - diagram:
            if all(tuple(a - int(j == i) for j, a in enumerate(v)) in diagram
                   for i in range(x) if v[i] > 0):
                result.add(diagram | {v})
    return result


def active(d, x):
    return all(any(v[i] > 0 for v in d) for i in range(x))


def image(d, permutation):
    return frozenset(tuple(v[i] for i in permutation) for v in d)


def check_orbits(family, m, x):
    remaining = set(family)
    number_orbits = 0
    orbit_sum = 0
    while remaining:
        d = min(remaining, key=lambda s: tuple(sorted(s)))
        q = frozenset(v for v in d if sum(v) == 2)
        h = len(q)
        e = sum(sum(a > 0 for a in v) == 2 for v in q)
        assert active(q, x), ("degree-two support", m, x, d)
        assert x <= 2 * h <= 2 * m
        orbit, stabilizer, kernel, q_images = set(), 0, 0, set()
        kernel_restrictions = set()
        q_order = sorted(q)
        for permutation in permutations(range(x)):
            transformed = image(d, permutation)
            orbit.add(transformed)
            if transformed == d:
                stabilizer += 1
                q_action = tuple(tuple(v[i] for i in permutation) for v in q_order)
                q_images.add(q_action)
                if q_action == tuple(q_order):
                    kernel += 1
                    restriction = []
                    for v in q_order:
                        support = [i for i, a in enumerate(v) if a]
                        if len(support) == 1:
                            assert permutation[support[0]] == support[0]
                        else:
                            assert set(permutation[i] for i in support) == set(support)
                            restriction.append(permutation[support[0]] == support[1])
                    kernel_restrictions.add(tuple(restriction))
        assert orbit <= family
        assert len(orbit) * stabilizer == factorial(x)
        assert len(kernel_restrictions) == kernel
        assert 2 ** e % kernel == 0
        assert factorial(h) % len(q_images) == 0
        assert kernel * len(q_images) == stabilizer
        assert (2 ** e * factorial(h)) % stabilizer == 0
        denominator = 2 ** m * factorial(m)
        assert denominator % stabilizer == 0
        orbit_sum += denominator // stabilizer
        remaining -= orbit
        number_orbits += 1
    assert Fraction(2 ** m * factorial(m) * len(family), factorial(x)) == orbit_sum
    return number_orbits


def polynomial(row):
    """Ascending coefficient list in Q[t], assembled in falling-factorial basis."""
    answer = [Fraction(0)] * len(row)
    falling = [1]
    for x, count in enumerate(row):
        for j, coefficient in enumerate(falling):
            answer[j] += Fraction(count, factorial(x)) * coefficient
        next_falling = [0] * (len(falling) + 1)
        for j, coefficient in enumerate(falling):
            next_falling[j] -= x * coefficient
            next_falling[j + 1] += coefficient
        falling = next_falling
    return answer


def evaluate(p, r):
    return sum(coefficient * r ** j for j, coefficient in enumerate(p))


def main():
    counts = [[0] * (2 * m + 1) for m in range(MAX_M + 1)]
    a = {}
    orbit_counts = [0] * (MAX_M + 1)
    for x in range(2 * MAX_M + 1):
        standard = mu(x)
        diagrams = {standard}
        for m in range(MAX_M + 1):
            a[m, x] = len(diagrams)
            family = {diagram - standard for diagram in diagrams
                      if active(diagram - standard, x)}
            if x <= 2 * m:
                counts[m][x] = len(family)
                orbit_counts[m] += check_orbits(family, m, x)
            else:
                assert not family
            if m != MAX_M:
                diagrams = grow(diagrams, x)
    assert counts == SOURCE_C, ("source C-matrix fidelity", counts)
    print("FIDELITY: Section 2.5 C-matrix rows m=0,...,4 reproduced by exhaustive labelled enumeration")
    for m, row in enumerate(counts):
        p = polynomial(row)
        denominator = 2 ** m * factorial(m)
        scaled = [denominator * coefficient for coefficient in p]
        assert all(coefficient.denominator == 1 for coefficient in scaled)
        assert scaled[-1] == 1
        if m:
            assert scaled[0] == 0
        for r in range(2 * MAX_M + 1):
            transformed = sum(comb(r, x) * count for x, count in enumerate(row) if x <= r)
            assert transformed == a[m, r] == evaluate(p, r)
        assert row[-1] == factorial(2 * m) // denominator
        print(f"m={m}: c={row}; alpha(z ascending)={[int(c) for c in reversed(scaled)]}; orbits={orbit_counts[m]}")
    assert [2 * c for c in polynomial(counts[1])] == [0, 1, 1]
    assert [8 * c for c in polynomial(counts[2])] == [0, 6, -1, 2, 1]
    print("SETTLEMENT CHECK: all strict skew diagrams for 0<=m<=4, 0<=x<=2m; degree-two coverage, kernel injection, stabilizer divisibility, orbit cancellation and integral coefficients")
    print("TRANSFORM CHECK: independently enumerated a_(m+r+1,r) for 0<=m<=4, 0<=r<=8")
    print("No unbounded assertion is certified by these finite checks.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
