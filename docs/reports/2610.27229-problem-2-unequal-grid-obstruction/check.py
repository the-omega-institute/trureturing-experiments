#!/usr/bin/env python3
"""Exact finite checks for García Guirao, arXiv:2610.27229v1, Problem 2.

The finite checks test the proposed obstruction, not the full ω-limit theorem.
Only Python integers and fractions.Fraction enter the arithmetic.
"""
from fractions import Fraction as Q
from itertools import product
from math import factorial


class PL:
    """Continuous linear interpolation of rational knot data on [0,1]."""

    def __init__(self, knots, values):
        self.knots, self.values = tuple(knots), tuple(values)
        assert self.knots[0] == 0 and self.knots[-1] == 1
        assert len(self.knots) == len(self.values)
        assert all(a < b for a, b in zip(self.knots, self.knots[1:]))
        assert all(Q(0) <= v <= Q(1) for v in self.values)

    def __call__(self, t):
        assert isinstance(t, Q) and 0 <= t <= 1
        for a, b, u, v in zip(self.knots, self.knots[1:], self.values, self.values[1:]):
            if a <= t <= b:
                return u + (v - u) * (t - a) / (b - a)
        raise AssertionError("outside interpolation domain")

    def range(self):
        return min(self.values), max(self.values)


def iterate(f, n, z):
    for _ in range(n):
        z = f(z)
    return z


def cournot(f1, f2, z):
    x, y = z
    return f2(y), f1(x)


def source_fidelity():
    # Introduction: F^(2n)(x,y) = ((f2∘f1)^n(x), (f1∘f2)^n(y)).
    grid = (Q(0), Q(1, 2), Q(1))
    maps = [PL(grid, values) for values in product(grid, repeat=3)]
    count = 0
    for f1, f2 in product(maps, repeat=2):
        F = lambda z: cournot(f1, f2, z)
        phi = lambda t: f2(f1(t))
        chi = lambda t: f1(f2(t))
        for x, y in product(grid, repeat=2):
            for n in range(7):
                assert iterate(F, 2 * n, (x, y)) == (iterate(phi, n, x), iterate(chi, n, y))
                count += 1
    print(f"FIDELITY even-iterate identity: {count} exact equalities; 729 map pairs, 9 points, n=0..6")

    # Section 10, Example 'A robust lattice': the displayed three-piece map g.
    g = PL((Q(0), Q(1, 2), Q(3, 4), Q(1)), (Q(1, 4), Q(1, 2), Q(1), Q(1, 2)))
    count = 0
    for t in (Q(0), Q(1, 16), Q(1, 4), Q(7, 16)):
        for n in range(21):
            expected = Q(1, 2) + (t - Q(1, 2)) / (2 ** n)
            assert iterate(g, n, t) == expected
            count += 1
    assert g(Q(1, 2)) == Q(1, 2)
    assert (min(g(t) for t in (Q(1, 2), Q(3, 4), Q(1))),
            max(g(t) for t in (Q(1, 2), Q(3, 4), Q(1)))) == (Q(1, 2), Q(1))
    print(f"FIDELITY robust-lattice map: {count} contraction equalities; g([1/2,1])=[1/2,1] exactly")


def in_H(z):
    x, y = z
    return 0 <= x <= 1 and 0 <= y <= 1 and (x in (0, 1) or y == Q(1, 2))


def in_image_H(f1, f2, z):
    # This is the exact union of the three segment images, not a sample of F(H).
    x, y = z
    a, b = f1.range()
    c, d = f2.range()
    return (x == f2(Q(1, 2)) and a <= y <= b) or (c <= x <= d and y in (f1(Q(0)), f1(Q(1))))


def obstruction_checks():
    totals = []
    for size in (3, 4):
        grid = tuple(Q(i, size - 1) for i in range(size))
        maps = [PL(grid, values) for values in product(grid, repeat=size)]
        # Three distinct y-values suffice to avoid the two horizontal image levels.
        candidates = (Q(0), Q(1, 2), Q(1))
        count = 0
        for f1, f2 in product(maps, repeat=2):
            y = next(t for t in candidates if t not in (f1(Q(0)), f1(Q(1))))
            x = Q(1) if f2(Q(1, 2)) == 0 else Q(0)
            z = x, y
            assert in_H(z)
            assert not in_image_H(f1, f2, z)
            count += 1
        totals.append(count)
        print(f"OBSTRUCTION H not contained in F(H): all {count} continuous PL map pairs on the {size}-point rational knot/value grid")
    assert totals == [729, 65536]

    # For a strongly invariant grid (P×J)∪(K×Q), finite transversal images
    # must satisfy f2(Q) ⊇ P and f1(P) ⊇ Q. Exhaust all finite assignments
    # between labelled P and Q for small cardinalities. These necessary
    # conditions alone are not asserted to construct an ω-limit realization.
    counts = []
    assignments = 0
    for m in range(1, 6):
        row = []
        for n in range(1, 6):
            P, R = set(range(m)), set(range(n))
            onto_qp = sum(set(v) == P for v in product(range(m), repeat=n))
            onto_pq = sum(set(v) == R for v in product(range(n), repeat=m))
            assignments += m ** n + n ** m
            possible = onto_qp * onto_pq
            assert (possible > 0) == (m == n)
            assert possible == (factorial(m) ** 2 if m == n else 0)
            row.append(possible)
        counts.append(row)
    print(f"GRID necessary two-surjection constraints: all {assignments} assignments, 1<=|P|,|Q|<=5")
    for m, row in enumerate(counts, 1):
        print(f"  |P|={m}, counts for |Q|=1..5: " + " ".join(map(str, row)))
    print("SCOPE finite exact checks only; no dense-orbit or general omega-limit theorem is certified here")


def main():
    source_fidelity()
    obstruction_checks()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
