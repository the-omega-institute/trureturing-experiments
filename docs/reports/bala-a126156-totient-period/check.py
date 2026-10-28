#!/usr/bin/env python3
"""Exact checks for Peter Bala's OEIS A126156 eventual totient-period conjecture.

Only Python 3 standard-library integer and rational arithmetic is used.
Source: https://oeis.org/A126156, revision 110; its integer triangle is
https://oeis.org/A126155, revision 26. OEIS data: CC BY-SA 4.0.
The finite checks corroborate the preregistered argument; they are not an
all-index proof or a Lean proof. Run: python3 check.py
"""
from fractions import Fraction
from math import comb, factorial, gcd

LIMIT = 200
SOURCE_TERMS = (
    1, 1, 7, 139, 5473, 357721, 34988647, 4784061619,
    871335013633, 203906055033841, 59618325600871687,
    21297483077038703899, 9127322584507530151393,
    4621897483978366951337161, 2730069675607609356178641127,
    1860452328661957054823447670979, 1448802510679254790311316267306753,
)
SOURCE_MOD55 = (
    1, 1, 7, 29, 28, 1, 12, 29, 18, 6, 12, 34, 18, 51,
    17, 34, 23, 51, 7, 39, 23, 1, 7, 29, 28, 1, 12, 29,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def triangle_rows(limit):
    """Literal zero-based A126155 PARI recurrence, including reflection."""
    rows = [[1]]
    for n in range(1, limit + 1):
        previous = rows[-1]
        first = sum(previous)
        lower = [first, 5 * first]
        for j in range(n - 1):
            lower.append(2 * lower[j + 1] - lower[j] - 8 * previous[j])
        rows.append(lower + lower[-2::-1])
    return rows


def tangent_recurrence(limit):
    """Tangent derivative recurrence and F' = F*tan(sqrt(2)*x)/sqrt(2)."""
    tangent = [1]
    for r in range(1, limit):
        tangent.append(sum(comb(2 * r, 2 * j + 1) * tangent[j]
                           * tangent[r - 1 - j] for j in range(r)))
    values = [1]
    for n in range(limit):
        values.append(sum(comb(2 * n + 1, 2 * m + 1) * (2 ** m)
                          * tangent[m] * values[n - m] for m in range(n + 1)))
    return values, tangent


def defining_egf(limit):
    """Recover F from F^2*cos(sqrt(2)*x)=1, F(0)=1, over Q[[x^2]].

    First invert cos, then solve the square root coefficient by coefficient.
    No symbolic trigonometric simplification or floating point is used.
    """
    cosine = [Fraction((-2) ** n, factorial(2 * n)) for n in range(limit + 1)]
    secant = [Fraction(1)]
    root = [Fraction(1)]
    values = [1]
    for n in range(1, limit + 1):
        secant.append(-sum(cosine[j] * secant[n - j] for j in range(1, n + 1)))
        root.append((secant[n] - sum(root[j] * root[n - j]
                                    for j in range(1, n))) / 2)
        value = root[n] * factorial(2 * n)
        require(value.denominator == 1, f"nonintegral e.g.f. coefficient at {n}")
        values.append(value.numerator)
    return values


def main():
    rows = triangle_rows(LIMIT)
    a = [row[0] for row in rows]
    tangent_values, tangent = tangent_recurrence(LIMIT)
    egf_values = defining_egf(LIMIT)
    require(a == tangent_values == egf_values, "three exact constructions disagree")
    print(f"Exact constructions: A126155 triangle, tangent recurrence, defining e.g.f.; agree at n=0..{LIMIT}")

    require(tuple(a[:len(SOURCE_TERMS)]) == SOURCE_TERMS, "OEIS term fidelity failed")
    require(rows[4] == [5473, 27365, 48145, 63365, 69025, 63365, 48145, 27365, 5473],
            "OEIS A126155 displayed row failed")
    require(tuple(v % 55 for v in a[:len(SOURCE_MOD55)]) == SOURCE_MOD55,
            "Bala's displayed mod-55 example failed")
    require(all((a[n + 20] - a[n]) % 55 == 0 for n in range(1, LIMIT - 19)),
            "mod-55 period-20 small-case check failed")
    print("Fidelity: all 17 displayed A126156 terms and the displayed A126155 row n=4 agree")
    print("Fidelity: all 28 displayed mod-55 residues agree; shift 20 agrees for n=1..180")

    require(tangent[0] == 1 and tangent[1] == 2, "tangent initial values failed")
    require(all((2 ** m * tangent[m]) % 4 == 0 for m in range(1, LIMIT)),
            "tangent divisibility failed")
    for n, row in enumerate(rows):
        require(len(row) == 2 * n + 1 and row == row[::-1], f"row shape failed at {n}")
        require(all(v % 4 == a[n] % 4 for v in row), f"row invariant failed at {n}")
    require(all((a[n + 1] - (2 * n + 1) * a[n]) % 4 == 0 for n in range(LIMIT)),
            "one-step congruence failed")
    require(all(v % 2 == 1 for v in a), "oddness failed")
    require(all((a[n + 2] + a[n]) % 4 == 0 for n in range(LIMIT - 1)),
            "two-step negative congruence failed")
    require(all((a[n + 2] - a[n]) % 4 == 2 for n in range(LIMIT - 1)),
            "two-step obstruction failed")
    require(all(v % 4 == (1 if (n // 2) % 2 == 0 else 3) for n, v in enumerate(a)),
            "four-cycle failed")
    print(f"Mod 4: every entry of every row n=0..{LIMIT} equals its leading entry")
    print("Mod 4: one-step congruence n=0..199; oddness n=0..200; two-step obstruction n=0..198")
    print("Mod-4 residues n=0..15:", [v % 4 for v in a[:16]])
    phi4 = sum(gcd(j, 4) == 1 for j in range(1, 5))
    periods = [p for p in range(1, phi4 + 1) if phi4 % p == 0]
    require(phi4 == 2 and periods == [1, 2], "totient/divisor calculation failed")
    # For every tested onset and every allowed p, locate an explicit failed shift.
    for onset in range(LIMIT - 2):
        for p in periods:
            witnesses = [n for n in range(onset, LIMIT - p + 1)
                         if (a[n + p] - a[n]) % 4 != 0]
            require(bool(witnesses), f"no obstruction for onset={onset}, p={p}")
    print("phi(4)=2; allowed periods [1,2]; failed shifts found for each onset N=0..197")
    print(f"Concrete positive indices: a(2)={a[2]}, a(4)={a[4]}, difference={a[4] - a[2]}, residue={(a[4] - a[2]) % 4}")
    print("Scope: finite exact evidence; arbitrary onset excluded by the preregistered all-index argument")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
