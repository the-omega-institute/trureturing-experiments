#!/usr/bin/env python3
"""Exact finite checks for Fish, arXiv:1702.02544v2, Conjecture 1.

These computations test a proposed ordinary proof; they are not a Lean proof
or an exhaustive prior-art certificate. Only Python's standard library is used.
"""
from fractions import Fraction
from itertools import product
from math import comb


def subsets(n):
    for mask in range(1 << n):
        yield {i for i in range(n) if (mask >> i) & 1}


def source_fidelity():
    # Section 3, proof of Corollary 1.2: R(E) = (E-E) + NZ.
    cases = 0
    dense_differences = 0
    for n in range(1, 10):
        for e in subsets(n):
            differences = {(a - b) % n for a in e for b in e}
            for t in range(-2 * n, 2 * n + 1):
                intersection = e & {(b + t) % n for b in e}
                measure = Fraction(len(intersection), n)
                assert (measure > 0) == (t % n in differences)
                cases += 1
            if 2 * len(e) > n:
                assert differences == set(range(n))
                # High-density specialization of Corollary 1.2, d=1.
                assert {a * b % n for a in differences for b in differences} == set(range(n))
                dense_differences += 1
    print(f"SOURCE_FIDELITY: return-time identity {cases} cases; dense cyclic sets {dense_differences}")


def boolean_balls():
    cases = 0
    for n in range(0, 9):
        universe = range(1 << n)
        balls = [{x for x in universe if x.bit_count() <= r} for r in range(n + 1)]
        for r, e in enumerate(balls):
            d = {x ^ y for x in e for y in e}
            assert d == balls[min(2 * r, n)]
            p = {x & y for x in d for y in d}
            assert p == d
            # Membership by four independent variables, not a generated ring.
            if n <= 5:
                p4 = {(a ^ b) & (c ^ f) for a, b, c, f in product(e, repeat=4)}
                assert p4 == p
            cases += 1
    print(f"BOOLEAN_BALLS: D(B_r)=B_min(2r,n), D*D=D in {cases} cases, n<=8")


def additive_subspaces(n):
    spaces = {frozenset({0})}
    frontier = set(spaces)
    while frontier:
        new = set()
        for s in frontier:
            for x in range(1 << n):
                if x not in s:
                    t = frozenset(set(s) | {x ^ y for y in s})
                    if t not in spaces:
                        new.add(t)
        spaces |= new
        frontier = new
    return spaces


def all_boolean_subrings():
    counts = []
    for n in range(0, 5):
        rings = [s for s in additive_subspaces(n) if all(x & y in s for x in s for y in s)]
        counts.append(len(rings))
        for s in rings:
            join = 0
            for x in s:
                join |= x
            assert join in s
            assert all(x & join == x for x in s)
        for bound in range(n + 1):
            contained = [s for s in rings if all(x.bit_count() <= bound for x in s)]
            largest = max(map(len, contained))
            assert largest == 1 << bound
            assert min(Fraction(1 << n, len(s)) for s in contained) == 1 << (n - bound)
    assert counts == [1, 2, 5, 15, 52]
    print(f"ALL_BOOLEAN_SUBRINGS: counts n=0..4 {counts}; exact support-size optimum checked")


def central_binomial_and_density():
    for t in range(513):
        p = Fraction(comb(2 * t, t), 4 ** t)
        q = Fraction(comb(2 * t + 2, t + 1), 4 ** (t + 1))
        assert q == p * Fraction(2 * t + 1, 2 * t + 2)
        assert p * p * (t + 1) <= 1
        assert 4 * (t + 1) ** 3 - (2 * t + 1) ** 2 * (t + 2) == 3 * t + 2
    print("CENTRAL_BINOMIAL: recurrence, square bound, induction polynomial for t=0..512")
    for m in range(2, 9):
        M = m ** 4
        n, r = 2 * M, M - m
        coefficient, count = 1, 0
        for j in range(r + 1):
            count += coefficient
            coefficient = coefficient * (n - j) // (j + 1)
        density = Fraction(count, 1 << n)
        p = Fraction(comb(2 * M, M), 4 ** M)
        symmetry = Fraction(1, 2) - p / 2 - sum(
            (Fraction(comb(2 * M, M - j), 4 ** M) for j in range(1, m)), Fraction(0)
        )
        assert density == symmetry
        assert p * p * (M + 1) <= 1
        assert p * m * m <= 1
        assert density >= Fraction(1, 2) - m * p >= Fraction(1, 2) - Fraction(1, m)
        assert 2 * r == n - 2 * m
        assert Fraction(1 << n, 1 << (n - 2 * m)) == 1 << (2 * m)
        print(f"DENSITY_FAMILY: m={m}, n={n}, r={r}, density>=1/2-1/{m}, optimal index={1 << (2*m)}")


def endpoint_family():
    for d in range(1, 5):
        # Encode (a,v,c) by a + 2*v + 2^(d+1)*c; addition is xor.
        size = 1 << (d + 2)
        z = 1 << (d + 1)
        mask = (1 << d) - 1
        def multiply(x, y):
            a, v, c = x & 1, (x >> 1) & mask, (x >> (d + 1)) & 1
            b, w, h = y & 1, (y >> 1) & mask, (y >> (d + 1)) & 1
            return (a & b) | (((a * w) ^ (b * v)) << 1) | (
                ((a * h) ^ (b * c) ^ ((v & w).bit_count() & 1)) << (d + 1)
            )
        table = [[multiply(x, y) for y in range(size)] for x in range(size)]
        for x, y, t in product(range(size), repeat=3):
            assert table[table[x][y]][t] == table[x][table[y][t]]
            assert table[x][y ^ t] == table[x][y] ^ table[x][t]
        assert all(table[x][1] == x and table[x][0] == 0 for x in range(size))
        assert all(table[x][y] == table[y][x] for x in range(size) for y in range(size))
        ideal = {x for x in range(size) if not x & 1}
        units = {x for x in range(size) if any(table[x][y] == 1 for y in range(size))}
        assert units == {1 ^ x for x in ideal}
        assert Fraction(len(units), size) == Fraction(1, 2)
        assert {x ^ y for x in units for y in units} == ideal
        square = {table[x][y] for x in ideal for y in ideal}
        assert square == {0, z}
        assert all(table[table[x][y]][t] == 0 for x, y, t in product(ideal, repeat=3))
        assert all(table[1 ^ x][1 ^ x ^ table[x][x]] == 1 for x in ideal)
        # Enumerate every subring of the actual product set.
        candidates = [set(), {0}, {z}, {0, z}]
        subrings = [s for s in candidates if 0 in s
                    and all(x ^ y in s and table[x][y] in s for x in s for y in s)]
        assert max(map(len, subrings)) == 2
        assert Fraction(size, 2) == 1 << (d + 1)
        print(f"ENDPOINT: d={d}, |A|={size}, density=1/2, |m*m|=2, minimum index={size//2}; full ring-law checks")


def main():
    assert __debug__, "Run without Python -O"
    source_fidelity()
    boolean_balls()
    all_boolean_subrings()
    central_binomial_and_density()
    endpoint_family()
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
