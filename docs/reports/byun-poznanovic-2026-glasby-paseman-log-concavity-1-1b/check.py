#!/usr/bin/env python3
"""Exact bounded checks for arXiv:2604.14639v1, Conjecture 1.1(b).

Only integers and Fraction are used. No bounded check establishes eventual
log-concavity for all positive real weights and all positive integer powers.
"""
from fractions import Fraction
from functools import lru_cache
from math import comb


@lru_cache(maxsize=None)
def denominator_rows(u, v, power, stop):
    """Return v**(r*power) times the literal denominator, for 0 <= r <= stop."""
    up = [u ** (i * power) for i in range(stop + 1)]
    vp = [v ** (i * power) for i in range(stop + 1)]
    return tuple(sum(comb(r, i) ** power * up[i] * vp[r-i]
                     for i in range(r + 1)) for r in range(stop + 1))


def scaled_numerators(m, u, v, power):
    # Both numerator and denominator are multiplied by the same v**(r*power).
    acc = 0
    result = []
    for r in range(m + 1):
        acc = v**power * acc + comb(m, r)**power * u**(r*power)
        result.append(acc)
    return result


def bad_indices(m, u, v, power, rows):
    ns = scaled_numerators(m, u, v, power)
    return [r for r in range(1, m)
            if ns[r]**2 * rows[r-1] * rows[r+1]
            < ns[r-1] * ns[r+1] * rows[r]**2]


def literal(m, r, a, power):
    return (sum((comb(m, i)*a**i)**power for i in range(r+1)) /
            sum((comb(r, i)*a**i)**power for i in range(r+1)))


def main():
    weights = (Fraction(1, 7), Fraction(1, 2), Fraction(1),
               Fraction(3, 2), Fraction(3), Fraction(10))
    comparisons = 0
    for a in weights:
        u, v = a.numerator, a.denominator
        for power in range(1, 5):
            rows = denominator_rows(u, v, power, 8)
            for m in range(9):
                ns = scaled_numerators(m, u, v, power)
                for r in range(m+1):
                    assert Fraction(ns[r], rows[r]) == literal(m, r, a, power)
                    comparisons += 1
    print(f"Literal Fraction/scaled-integer agreement: {comparisons} entries")

    # Source fidelity: equation (1.5), Theorem 1.2(a),(b), and Table 2.
    for r in range(101):
        assert sum(comb(r, i)**2 for i in range(r+1)) == comb(2*r, r)
    for m in range(2, 101):
        rows = denominator_rows(1, 1, 2, 100)
        assert not bad_indices(m, 1, 1, 2, rows)
        ns = scaled_numerators(m, 1, 1, 2)
        peak = (m+2)//3
        assert all(ns[peak]*rows[r] > ns[r]*rows[peak]
                   for r in range(m+1) if r != peak)
    print("Source eq. (1.5): central-binomial identity r=0..100")
    print("Source Theorem 1.2(a),(b): log-concavity and unique peak m=2..100")

    # Literal transcription of Table 2: rows l=1..20, columns a=1..10.
    table = [[0]*10 for _ in range(5)] + [
        [5,0,0,0,0,0,0,0,0,0], [5,0,0,0,0,0,0,0,0,0],
        [9,0,0,0,0,0,0,0,0,0], [12,7,10,0,0,0,0,0,0,0],
        [13,7,10,0,0,0,0,0,0,0], [16,10,15,13,15,0,0,0,0,0],
        [19,12,15,18,21,24,28,0,0,0], [20,14,20,19,21,25,28,32,35,39],
        [24,16,20,24,28,25,29,32,36,39],
        [27,17,25,25,28,33,37,42,46,40],
        [30,20,25,25,29,33,38,42,47,51],
        [31,21,30,31,35,41,46,42,47,52],
        [35,25,30,31,36,41,46,52,57,63],
        [39,25,35,37,42,41,47,52,58,63],
        [42,29,35,37,43,49,55,62,68,64],
    ]
    assert len(table) == 20
    cases = 0
    for power, expected in enumerate(table, 1):
        for a, last in enumerate(expected, 1):
            rows = denominator_rows(a, 1, power, 100)
            failures = [m for m in range(1, 101)
                        if bad_indices(m, a, 1, power, rows)]
            assert max(failures, default=0) == last, (power, a, failures, last)
            cases += 100
    print(f"Source Table 2: all 200 cells reproduced, {cases} complete sequences")

    # Proposed analytic route: all-index denominator bound and strict finite-head limit.
    bound_count = 0
    for a in weights:
        u, v = a.numerator, a.denominator
        for power in range(1, 13):
            ds = denominator_rows(u, v, power, 101)
            for r in range(1, 32):
                assert ds[r-1]*ds[r+1]*(r+2)**power >= ds[r]**2*(r+1)**power
                assert ds[r-1]*ds[r+1]*(r+1)**power > ds[r]**2*r**power
                bound_count += 1
    print(f"Denominator lower curvature bound and strict fixed-r limit: {bound_count} triples")

    # An additional nonintegral-weight check; finite exceptions are recorded, not denied.
    checked = 0
    last_failures = []
    for a in weights:
        u, v = a.numerator, a.denominator
        worst = 0
        for power in range(1, 13):
            ds = denominator_rows(u, v, power, 101)
            failures = [m for m in range(1, 101)
                        if bad_indices(m, u, v, power, ds)]
            last = max(failures, default=0)
            worst = max(worst, last)
            for m in range(1, 101):
                ns = scaled_numerators(m, u, v, power)
                assert all(ns[r]**2 >= ns[r-1]*ns[r+1] for r in range(1, m))
            checked += 100
        last_failures.append((str(a), worst))
    print(f"Additional grid: {checked} sequences, l=1..12, m=1..100")
    print(f"Largest observed exceptional m over l=1..12, by weight: {last_failures}")
    assert 18 in bad_indices(100, 1, 7, 10, denominator_rows(1, 7, 10, 101))
    assert literal(100, 18, Fraction(1, 7), 10)**2 < (
        literal(100, 17, Fraction(1, 7), 10) * literal(100, 19, Fraction(1, 7), 10))
    print("Numerator prefixes log-concave on the entire additional grid")
    print("Exact finite exception retained: a=1/7, l=10, m=100, r=18")
    print("A largest observed exceptional m of 100 is censored by the scan ceiling")
    print("Scope: bounded exact evidence only; the universal asymptotic claim is not tested by enumeration")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
