#!/usr/bin/env python3
"""Exact, bounded checks of Mathar's OEIS A084569 square interpretation.

Independent enumeration of all eight positive entries for n <= 4 is compared
with enumeration after eliminating the four side-middle entries for n <= 16.
This program is finite evidence, not a proof of the universally quantified claim.
No external packages, network calls, randomness, or floating point are used.
"""
from collections import Counter
from fractions import Fraction
from itertools import product

# OEIS A084569 revision 19, DATA; README v0.069, Example 7, same prefix.
SOURCE_PREFIX = (1, 3, 9, 21, 44, 82, 142, 230, 355)
# OEIS A084569 COMMENTS, read counterclockwise from a corner, magic sum 4.
SOURCE_SUM4 = (
    (1, 1, 2, 1, 1, 1, 2, 1),
    (1, 1, 2, 1, 1, 2, 1, 2),
    (1, 2, 1, 2, 1, 2, 1, 2),
)


def shift(v, k):
    return v[k:] + v[:k]


def square_orbit(v):
    return frozenset(shift(v, 2 * k) for k in range(4))


def is_square(v, n):
    return all(t >= 1 for t in v) and all(
        v[2 * i] + v[2 * i + 1] + v[(2 * i + 2) % 8] == n + 3
        for i in range(4)
    )


def full_enumeration(n):
    # A side sum is n+3 and its other two labels are >=1, so each label <=n+1.
    return {v for v in product(range(1, n + 2), repeat=8) if is_square(v, n)}


def corners(n):
    return {x for x in product(range(n + 1), repeat=4)
            if all(x[i] + x[(i + 1) % 4] <= n for i in range(4))}


def to_square(x, n):
    return tuple(t for i in range(4)
                 for t in (x[i] + 1, n - x[i] - x[(i + 1) % 4] + 1))


def source_sum(n):
    return sum(sum(sum(i + (-1) ** i for i in range(j + 1))
                   for j in range(k + 1)) for k in range(n + 1))


def source_formula(n):
    return Fraction((-1) ** n, 8) + Fraction(
        n**4 + 6*n**3 + 17*n**2 + 30*n + 21, 24)


def main():
    source_orbits = {square_orbit(v) for v in SOURCE_SUM4}
    assert all(is_square(v, 1) for v in SOURCE_SUM4)
    assert len(source_orbits) == 3
    assert source_orbits == {square_orbit(v) for v in full_enumeration(1)}
    print('FIDELITY: all three magic-sum-4 squares in OEIS COMMENTS reproduced')

    for n in range(5):
        full = full_enumeration(n)
        reduced = {to_square(x, n) for x in corners(n)}
        assert full == reduced, (n, 'eight-entry vs corner enumeration')
        assert all(tuple(v[2*i] - 1 for i in range(4)) in corners(n) for v in full)
    print('FIDELITY: eight-entry enumeration equals corner reconstruction for n=0..4')

    print('n magic_sum E H Q orbit_count stabilizer1 stabilizer2 stabilizer4')
    observed = []
    for n in range(17):
        xs = corners(n)
        fixes = tuple(sum(shift(x, k) == x for x in xs) for k in range(4))
        orbits = {frozenset(shift(x, k) for k in range(4)) for x in xs}
        square_orbits = {square_orbit(to_square(x, n)) for x in xs}
        assert len(square_orbits) == len(orbits)
        # Verify inverse reconstruction and equivariance for every enumerated tuple.
        for x in xs:
            v = to_square(x, n)
            assert is_square(v, n)
            assert tuple(v[2*i] - 1 for i in range(4)) == x
            for k in range(4):
                assert to_square(shift(x, k), n) == shift(v, 2*k)
        E, Q, H, Q3 = fixes
        assert E == sum((2*m + 1)*(n-m+1)**2 for m in range(n+1))
        assert E == Fraction((n+1)*(n+2)*(n*n+3*n+3), 6)
        assert H == Fraction((n+1)*(n+2), 2)
        assert Q == Q3 == n//2 + 1
        assert 4*len(orbits) == sum(fixes)
        assert len(orbits) == source_sum(n) == source_formula(n)
        sizes = Counter(4//len(o) for o in orbits)
        assert set(sizes) <= {1, 2, 4}
        assert sum(sizes.values()) == len(orbits)
        assert 4*sizes[1] + 2*sizes[2] + sizes[4] == E
        assert 2*sizes[2] + sizes[4] == H
        assert sizes[4] == Q
        observed.append(len(orbits))
        print(n, n+3, E, H, Q, len(orbits), sizes[1], sizes[2], sizes[4])

    assert tuple(observed[:len(SOURCE_PREFIX)]) == SOURCE_PREFIX
    print('FIDELITY: README Example 7 / OEIS DATA prefix reproduced exactly')
    for n in range(41):
        assert source_sum(n) == source_formula(n)
    print('SEQUENCE: source triple sum equals source quartic formula for n=0..40')

    # Distinguish the binding C4 quotient from an unintended D4 quotient.
    vs = {to_square(x, 3) for x in corners(3)}
    cyclic = {min(square_orbit(v)) for v in vs}
    dihedral = {min(square_orbit(v) | square_orbit(
        tuple(v[(-i) % 8] for i in range(8)))) for v in vs}
    assert len(cyclic) == 21 and len(dihedral) == 19
    # A shift by one perimeter position exchanges corners and middles and is not
    # a physical quarter-turn; it can even change the side sum.
    assert is_square(SOURCE_SUM4[2], 1)
    assert not is_square(shift(SOURCE_SUM4[2], 1), 1)
    print('CONVENTIONS: n=3 has 21 C4 orbits and 19 D4 orbits; odd shifts excluded')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
