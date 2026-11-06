#!/usr/bin/env python3
"""Exact checks for Devlin--Zeng's even-exponent density conjectures.

All counts use Python integers and Fraction; no floating point is used as evidence.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import permutations
from math import comb


def kaprekar_digits(base: int, digits: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    """The source's K_b on a four-digit string, retaining leading zeros."""
    down = tuple(sorted(digits, reverse=True))
    up = tuple(sorted(digits))
    d = sum(a * base ** (3 - i) for i, a in enumerate(down))
    a = sum(a * base ** (3 - i) for i, a in enumerate(up))
    z = d - a
    return tuple((z // base ** (3 - i)) % base for i in range(4))


def pair_map(base: int, pair: tuple[int, int]) -> tuple[int, int]:
    d, e = pair
    if d < e:
        d, e = e, d
    if d == 0 and e == 0:
        return (0, 0)
    if e == 0:
        return tuple(sorted((d - 1, base - d), reverse=True))
    if d == e or d + e == base:
        return tuple(sorted((abs(2 * d - (base - 1)), abs(2 * d - (base + 1))), reverse=True))
    return tuple(sorted((abs(2 * d - base), abs(2 * e - base)), reverse=True))


def pair_weight(base: int, pair: tuple[int, int]) -> int:
    """Number of ordered four-digit strings with the given difference pair."""
    d, e = pair
    if not (0 <= e <= d < base):
        return 0
    if d == 0:
        return base
    if e == d:
        return 6 * (base - d)
    if e == 0:
        return 4 * (base - d) * (3 * d - 1)
    return 24 * (base - d) * (d - e)


def reaches_fixed_pair(base: int, start: tuple[int, int]) -> bool:
    target = (3 * base // 5, base // 5)
    seen: set[tuple[int, int]] = set()
    p = start
    while p not in seen:
        if p == target:
            return True
        seen.add(p)
        p = pair_map(base, p)
    return p == target


def density(base: int) -> Fraction:
    total = 0
    for d in range(base):
        for e in range(d + 1):
            if reaches_fixed_pair(base, (d, e)):
                total += pair_weight(base, (d, e))
    return Fraction(total, base ** 4)


def direct_density(base: int) -> Fraction:
    """Small-base direct enumeration, used as a cross-check of pair weights."""
    hit = 0
    for ds in __import__("itertools").product(range(base), repeat=4):
        x = ds
        seen = set()
        while x not in seen:
            seen.add(x)
            y = kaprekar_digits(base, x)
            if y == x:
                break
            x = y
        if x == kaprekar_digits(base, x) and x != (0, 0, 0, 0):
            hit += 1
    return Fraction(hit, base ** 4)


def residue_family_count(base: int, residues: set[tuple[int, int]]) -> int:
    """Exact lower/upper witness count from the source's 10-residue construction."""
    M = base // 10
    # Four distinct quotient blocks q_i, each in [0,M), and a residue-pair
    # condition.  For each (d,e), r0,r1 are free and r2,r3 are determined.
    return len(residues) * 2400 * comb(M, 4)


def main() -> None:
    # Fidelity: the worked decimal trajectory in the source.
    expected = {
        (3, 2, 2, 3): (1, 0, 8, 9),
        (0, 8, 8, 9): (8, 9, 9, 1),
        (8, 9, 9, 1): (8, 0, 8, 2),
        (8, 0, 8, 2): (8, 5, 3, 2),
        (8, 5, 3, 2): (6, 1, 7, 4),
        (6, 1, 7, 4): (6, 1, 7, 4),
    }
    for x, y in expected.items():
        assert kaprekar_digits(10, x) == y, (x, kaprekar_digits(10, x), y)

    # The displayed first subtraction in the source has an arithmetic typo: 3322 - 2233 = 1089.
    # The following steps are exact, and the pair reduction below is independent of that typo.
    # Direct enumeration agrees with the pair reduction at bases where it is cheap.
    for b in (5, 10, 15, 20):
        assert direct_density(b) == density(b), (b, direct_density(b), density(b))

    # Fidelity to Devlin--Zeng Theorem (proportion): b=5*m*2^n with odd m>1.
    # For m=3,n=0, the stated value is (8+40)/(5*b^2)=16/375.
    assert density(15) == Fraction(16, 375), density(15)

    # Exact values requested by the preregistration orchestrator.
    values = {b: density(b) for b in (80, 320, 1280)}
    assert values[80] == Fraction(2355169, 2560000)
    assert values[320] == Fraction(246471427, 327680000)
    assert values[1280] == Fraction(38517137157, 41943040000)

    # Candidate's residue families.  For n=4r, m=2^n is 6 mod 10 and the
    # 76 ordered residue pairs in the proposed basin family give a 19/25 limit.
    basin_residues = {
        (d, e) for d in range(10) for e in range(10)
        if (d % 5 != 0 and e % 5 != 0)
        or ((d % 10 == 5 and e in {1, 2, 3, 7, 8, 9})
            or (e % 10 == 5 and d in {1, 2, 3, 7, 8, 9}))
    }
    assert len(basin_residues) == 76
    # For n=4r+2, the eight ordered residue pairs are trapped in a non-fixed cycle.
    obstructing = {(d, e) for d in {2, 8} for e in {3, 7}} | {
        (d, e) for d in {3, 7} for e in {2, 8}
    }
    assert len(obstructing) == 8
    b0 = 5 * 2 ** 4
    b1 = 5 * 2 ** 6
    b2 = 5 * 2 ** 8
    assert values[b0] > Fraction(3, 4)
    assert values[b1] < Fraction(15, 16)
    # The witness fractions converge to the claimed strict separation.
    for b, rs, target in ((b0, basin_residues, Fraction(19, 25)),
                          (b1, obstructing, Fraction(23, 25))):
        M = b // 10
        witness = Fraction(residue_family_count(b, rs), b ** 4)
        if rs is basin_residues:
            # The finite witness is below 19/25 and tends to 19/25 from below.
            assert witness < target
        else:
            # The eight classes are outside the basin, so C_b <= 1 - witness;
            # for b=320 this is already strictly below 15/16.
            assert Fraction(1, 1) - witness < Fraction(15, 16)
        assert M >= 4

    # Every pair in the 76-class family is attracted at b=80 and b=1280.
    for b, n in ((b0, 4), (b2, 8)):
        for d in range(b):
            for e in range(d + 1):
                if (d % 10, e % 10) in basin_residues:
                    assert reaches_fixed_pair(b, (d, e)), (b, d, e)

    # Explicit pair-map checks for the 8 obstruction residues at b=5*2^(4r+2).
    b = b1
    n = 6
    for r in obstructing:
        p = r
        for _ in range(n):
            p = pair_map(b, p)
        assert p == (3 * (b // 5), 2 * (b // 5)), (r, p)

    print("fidelity: decimal K_10 trajectory and Theorem (proportion) checked")
    for b, c in values.items():
        print(f"C_{b} = {c}  (against 3/4={Fraction(3,4)}, 15/16={Fraction(15,16)})")
    print("residue families: 76 attracted classes and 8 obstructing classes checked")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
