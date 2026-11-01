#!/usr/bin/env python3
"""Exact checks for the OEIS A126364/A025565 boundary refutation."""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path

import sympy as sp

A126364 = [
    10, 28, 80, 230, 664, 1922, 5574, 16188, 47064, 136946,
    398746, 1161634, 3385486, 9869934, 28781908, 83948652,
    244894048, 714493794, 2084792450, 6083620812, 17753709802,
    51812911858, 151218254846, 441351052720, 1288175317718,
    3759883395300, 10974394864558,
]
A025565 = [
    1, 2, 4, 10, 26, 70, 192, 534, 1500, 4246, 12092, 34606,
    99442, 286730, 829168, 2403834, 6984234, 20331558, 59287740,
    173149662, 506376222, 1482730098, 4346486256, 12754363650,
    37461564504, 110125172682, 323990062452,
]


def digit_counts(n: int) -> list[int]:
    """Counts by final digit of length-n words over 0..9 with adjacent gap <= 1."""
    if n < 1:
        raise ValueError("source sequence starts at n=1")
    row = [1] * 10
    for _ in range(1, n):
        row = [sum(row[max(0, d - 1): min(10, d + 2)]) for d in range(10)]
    return row


def bounded_walk_counts(r: int, upper: int | None) -> int:
    """Motzkin walks from 0, staying in [0,upper] (or [0,infinity])."""
    row = [1] + [0] * (upper if upper is not None else r + 1)
    for _ in range(r):
        nxt = [0] * len(row)
        for h, count in enumerate(row):
            if not count:
                continue
            for step in (-1, 0, 1):
                nh = h + step
                if 0 <= nh < len(row):
                    nxt[nh] += count
        row = nxt
    return sum(row)


def b_value(n: int) -> int:
    return 1 if n == 1 else 2 * bounded_walk_counts(n - 2, None)


def exact_gf_check() -> None:
    """Check the source generating function with exact SymPy coefficients."""
    x = sp.Symbol("x")
    gf = x * sp.sqrt((1 + x) / (1 - 3 * x))
    coeffs = [int(sp.expand(sp.series(gf, x, 0, 28).removeO()).coeff(x, n)) for n in range(1, 28)]
    assert coeffs == A025565
    # A Fraction sanity check keeps the algebra explicitly rational at the first nontrivial term.
    assert Fraction(coeffs[11], 1) == Fraction(34606, 1)


def main() -> None:
    exact_gf_check()
    a_values = [sum(digit_counts(n)) for n in range(1, 28)]
    b_values = [b_value(n) for n in range(1, 28)]
    assert a_values == A126364, (a_values, A126364)
    assert b_values == A025565, (b_values, A025565)
    print("fidelity: A126364 n=1..27 and A025565 n=1..27 reproduced exactly")

    for n in range(2, 201):
        a_n = sum(digit_counts(n))
        a_prev = sum(digit_counts(n - 1))
        b_n = b_value(n)
        u = bounded_walk_counts(n - 2, None)
        v = bounded_walk_counts(n - 2, 9)
        assert a_n - 3 * a_prev + b_n == 2 * (u - v)
        if n <= 11:
            assert a_n == 3 * a_prev - b_n
        else:
            assert a_n > 3 * a_prev - b_n
    print("defect identity and equality/failure threshold checked for every n=2..200")

    a11 = sum(digit_counts(11))
    a12 = sum(digit_counts(12))
    b12 = b_value(12)
    assert (a11, a12, b12) == (398746, 1161634, 34606)
    assert a12 - (3 * a11 - b12) == 2
    assert bounded_walk_counts(10, None) == 17303
    assert bounded_walk_counts(10, 9) == 17302
    print("n=12 witness: a(12)=1161634, 3*a(11)-b(12)=1161632, defect=2")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
