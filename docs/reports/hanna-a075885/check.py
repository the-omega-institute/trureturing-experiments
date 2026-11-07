#!/usr/bin/env python3
"""Exact diagnostics for the proposed A075885 limit proof; no limit inference."""

from collections import Counter, defaultdict
from fractions import Fraction
from math import prod


SOURCE_TERMS = (
    1, 2, 5, 10, 29, 46, 169, 239, 745, 1450, 4111, 5182,
    27157, 33164, 84001, 186496, 610065, 713474, 3061009,
    3526553, 13783421, 27380452, 63264389, 71240523, 444872761,
    620729126, 1400231613, 2615011102, 9094701085, 10008828958,
)
LIMIT = 1000


def prefixes(n):
    values = [1]
    for k in range(1, n + 1):
        values.append(values[-1] * (n // k))
    return values


def source_formula(n):
    # A separate implementation of the source's literal sum of products.
    return 1 + sum(prod(n // k for k in range(1, m + 1))
                   for m in range(1, n + 1))


def normalized_coefficients(coefficients):
    # These are coefficients of formal log(j), not approximate logarithms.
    # log(1) is exactly zero.
    return {j: c for j, c in coefficients.items() if j != 1 and c}


def main():
    for n, expected in enumerate(SOURCE_TERMS):
        assert source_formula(n) == expected, ("source term", n)
        assert sum(prefixes(n)) == expected, ("prefix recurrence", n)
    assert source_formula(5) == 46
    print("FIDELITY: all 30 source terms n=0..29 and the a(5)=46 example agree")

    for n in range(1, LIMIT + 1):
        ps = prefixes(n)
        p = ps[-1]
        a = sum(ps)
        h = n // 2
        c = n - h
        assert p > 0 and a > 0
        assert all(x == p for x in ps[h:]), ("plateau", n)
        for m in range(h):
            assert ps[m] * 2 ** (h - m) <= p, ("geometric bound", n, m)
        early = Fraction(sum(ps[:h]), p)
        assert 0 <= early <= 1 - Fraction(1, 2 ** h) < 1
        assert Fraction(a, p) == c + 1 + early
        assert (c + 1) * p <= a < (c + 2) * p
        assert p <= a <= (n + 1) * p

        counts = Counter(n // k for k in range(1, n + 1))
        assert all(counts[j] == n // j - n // (j + 1)
                   for j in range(1, n + 1)), ("floor count", n)
        grouped = prod(j ** (n // j - n // (j + 1))
                       for j in range(1, n + 1))
        assert p == grouped, ("grouped product", n)
        oudra = prod((Fraction(j, j - 1) ** (n // j)
                      for j in range(2, n + 1)), start=Fraction(1))
        assert oudra == p, ("Oudra rational product", n)

        lhs = defaultdict(Fraction)
        rhs = defaultdict(Fraction)
        for j in range(2, n + 1):
            lhs[j] += Fraction(1, j)
            lhs[j - 1] -= Fraction(1, j)
        for j in range(1, n):
            rhs[j] += Fraction(1, j * (j + 1))
        rhs[n] += Fraction(1, n)
        assert normalized_coefficients(lhs) == normalized_coefficients(rhs), (
            "summation by parts", n)

        for j in range(2, n + 1):
            error = Fraction(1, j) - Fraction(n // j, n)
            assert 0 <= error < Fraction(1, n), ("floor error", n, j)
        tail = sum((Fraction(1, j * (j - 1)) for j in range(n + 1, 2 * n + 1)),
                   start=Fraction(0))
        assert tail == Fraction(1, n) - Fraction(1, 2 * n)

    print(f"SETTLEMENT DIAGNOSTICS: n=1..{LIMIT}; exact integer/rational checks")
    print("PASS: terminal plateau, geometric early-mass bound, both sum/product sandwiches")
    print("PASS: floor multiplicities, grouped floor product, Oudra rational product identity")
    print("PASS: formal logarithm summation-by-parts coefficients, floor error, telescoping tail")
    print("The infinite limit and analytic logarithm inequalities are not inferred from these checks.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
