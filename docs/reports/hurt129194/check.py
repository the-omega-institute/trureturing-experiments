#!/usr/bin/env python3
"""Exact fidelity and finite checks for Hurt's OEIS A129194 conjecture."""
from fractions import Fraction
from math import gcd

# Data field of OEIS A129194 revision 94, offset 0,3 (CC BY-SA 4.0).
# Source: https://oeis.org/A129194 ; no OEIS program source is used.
DISPLAYED = (
    0, 1, 2, 9, 8, 25, 18, 49, 32, 81, 50, 121, 72, 169, 98, 225,
    128, 289, 162, 361, 200, 441, 242, 529, 288, 625, 338, 729, 392,
    841, 450, 961, 512, 1089, 578, 1225, 648, 1369, 722, 1521, 800,
    1681, 882, 1849, 968, 2025, 1058, 2209, 1152, 2401, 1250, 2601,
    1352,
)


def source_a(n: int) -> Fraction:
    """Literal NAME formula, interpreted in the rational numbers."""
    if n < 0:
        raise ValueError("OEIS indices start at zero")
    return Fraction(n, 2) ** 2 * (3 - (-1) ** n)


def rational_expression(n: int) -> Fraction:
    if n < 1:
        raise ValueError("the source rational expression requires n >= 1")
    return Fraction(4, n) - Fraction(2, n * n)


def check_index(n: int) -> None:
    q = rational_expression(n)
    a = source_a(n)
    assert a.denominator == 1, (n, a)
    assert q.denominator == a, (n, q, a)
    numerator, denominator = 2 * (2 * n - 1), n * n
    assert q == Fraction(numerator, denominator), (n, q)
    assert gcd(2 * n - 1, n * n) == 1, n
    assert gcd(numerator, denominator) == gcd(2, n * n), n
    assert q.denominator == (n * n) // gcd(2, n * n), n
    assert q.denominator == (n * n if n % 2 else (n * n) // 2), n


def main() -> None:
    for n, expected in enumerate(DISPLAYED):
        assert source_a(n) == expected, (n, source_a(n), expected)
    print(f"FIDELITY: all {len(DISPLAYED)} OEIS displayed terms n=0..52 match the literal formula")
    for n in range(6, 10001):
        assert source_a(n) == 3 * source_a(n - 2) - 3 * source_a(n - 4) + source_a(n - 6), n
    print("FIDELITY: source recurrence a(n)=3a(n-2)-3a(n-4)+a(n-6) matches for n=6..10000")
    for n in range(1, 10001):
        check_index(n)
    print("SETTLEMENT CHECK: denominator identity, fraction rewrite, gcd mechanism and parity formula match for every n=1..10000")
    for n in (10001, 10002, 10**30, 10**30 + 1):
        check_index(n)
    print("LARGE EXACT CHECKS: n=10001,10002,10^30,10^30+1 all match")
    for n in range(1, 7):
        print(f"EXAMPLE: n={n}, q={rational_expression(n)}, denominator={rational_expression(n).denominator}, a={source_a(n)}")
    try:
        rational_expression(0)
    except ValueError:
        print("BOUNDARY: a(0)=0; rational expression at n=0 is excluded as undefined")
    else:
        raise AssertionError("n=0 must be excluded")
    print("Finite exact checks support the proposed route; they do not prove the universal conjecture.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
