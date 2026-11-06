#!/usr/bin/env python3
"""Exact checks for the literal array count in OEIS A183624.

No recurrence or closed form is used to generate the array counts.
Exhaustion visits every labelled array, including inadmissible arrays.
The separate row-state dynamic program checks each four-entry block directly.
Finite checks support a preregistration; they do not prove the unbounded claim.
"""
from fractions import Fraction
from itertools import product


def count_exhaustively(m: int, rows: int) -> int:
    assert m >= 0 and rows >= 1
    return sum(
        all(
            cells[2*i] + cells[2*i+1] + cells[2*i+2] + cells[2*i+3] == 2*m
            for i in range(rows - 1)
        )
        for cells in product(range(m + 1), repeat=2*rows)
    )


def row_state_counts(m: int, max_rows: int) -> list[int]:
    states = list(product(range(m + 1), repeat=2))
    counts = {state: 1 for state in states}
    totals = [sum(counts.values())]
    for _ in range(1, max_rows):
        counts = {
            (u, v): sum(
                count for (s, t), count in counts.items()
                if s + t + u + v == 2*m
            )
            for u, v in states
        }
        totals.append(sum(counts.values()))
    return totals


def proposed_count(m: int, rows: int) -> int:
    return (m + 1)**rows + 2*sum(j**rows for j in range(1, m + 1))


def multiply(a: list[Fraction], b: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0) for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x*y
    return out


def main() -> None:
    checks = 0
    max_ambient = 0
    for m in range(5):
        dp = row_state_counts(m, 6)
        for rows in range(1, 7):
            ambient = (m + 1)**(2*rows)
            if ambient > 600_000:
                continue
            actual = count_exhaustively(m, rows)
            expected = proposed_count(m, rows)
            assert actual == dp[rows - 1] == expected, (m, rows, actual, expected)
            checks += 1
            max_ambient = max(max_ambient, ambient)
            print(f"EXHAUSTIVE m={m} rows={rows} arrays={ambient} count={actual} OK")
    assert checks == 26
    print(f"EXHAUSTIVE_SUMMARY checks={checks} max_arrays={max_ambient}")

    # Fidelity to the source's first five terms, via literal arrays.
    published = [19, 45, 115, 309, 859]
    measured = [count_exhaustively(2, n + 1) for n in range(1, 6)]
    assert measured == published, (measured, published)
    print(f"SOURCE_FIDELITY n=1..5 counts={measured} OK")

    # A second, exact enumerator checks n=1..40 without enumerating all arrays.
    dp = row_state_counts(2, 41)
    coefficients = [Fraction(0)] + [Fraction(dp[n]) for n in range(1, 41)]
    for n in range(1, 41):
        assert coefficients[n] == 2 + 2**(n + 2) + 3**(n + 1), n
    denominator = [Fraction(1)]
    for j in (1, 2, 3):
        denominator = multiply(denominator, [Fraction(1), Fraction(-j)])
    assert denominator == [1, -6, 11, -6]
    product_coefficients = multiply(coefficients, denominator)
    numerator = [Fraction(0), Fraction(19), Fraction(-69), Fraction(54)]
    for n in range(41):
        expected = numerator[n] if n < len(numerator) else Fraction(0)
        assert product_coefficients[n] == expected, (n, product_coefficients[n])
    print("NATIVE_ROW_DP closed_form_and_GF_coefficients n=0..40 OK")

    # Exact polynomial verification of the proposed geometric-series combination.
    factors = [[Fraction(1), Fraction(-j)] for j in (1, 2, 3)]
    partial_numerator = [Fraction(0) for _ in range(3)]
    for omitted, weight in enumerate((2, 8, 9)):
        term = [Fraction(1)]
        for i, factor in enumerate(factors):
            if i != omitted:
                term = multiply(term, factor)
        for i, c in enumerate(term):
            partial_numerator[i] += weight*c
    assert partial_numerator == [19, -69, 54]
    print("EXACT_POLYNOMIAL numerator_over_common_denominator=[19,-69,54] OK")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
