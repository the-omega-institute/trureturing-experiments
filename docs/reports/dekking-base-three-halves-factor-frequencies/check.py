#!/usr/bin/env python3
"""Exact finite checks for Dekking's base-3/2 factor-frequency conjecture.

All arithmetic used as evidence is integer or Fraction arithmetic.  The checks
are finite sanity checks for the definitions and the recurrence; they do not
claim to prove the limiting frequency statement.
"""

from fractions import Fraction


def sesquinary_digits(n: int) -> list[int]:
    """Return the canonical base-3/2 digits, least significant first."""
    if n < 0:
        raise ValueError("n must be nonnegative")
    if n == 0:
        return [0]
    out: list[int] = []
    while n:
        digit = n % 3
        out.append(digit)
        n = (2 * n - 2 * digit) // 3
    return out


def sesquinary_value(digits: list[int]) -> Fraction:
    value = Fraction(0)
    power = Fraction(1)
    for digit in digits:
        value += digit * power
        power *= Fraction(3, 2)
    return value


def x(n: int) -> int:
    return sum(sesquinary_digits(n)) % 2


def prefix(length: int) -> list[int]:
    return [x(n) for n in range(length)]


def count_factor(word: list[int], factor: tuple[int, ...]) -> int:
    m = len(factor)
    return sum(tuple(word[i : i + m]) == factor for i in range(len(word) - m + 1))


def main() -> None:
    # The paper's table and displayed prefix (A357448) are reproduced exactly.
    expected_digits = {
        0: [0],
        1: [1],
        2: [2],
        3: [0, 2],
        4: [1, 2],
        5: [2, 2],
        6: [0, 1, 2],
        7: [1, 1, 2],
        8: [2, 1, 2],
        9: [0, 0, 1, 2],
        10: [1, 0, 1, 2],
    }
    for n, digits in expected_digits.items():
        assert sesquinary_digits(n) == digits, (n, sesquinary_digits(n), digits)
        assert sesquinary_value(digits) == n, (n, sesquinary_value(digits))

    displayed = [0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0,
                 1, 0, 1, 0, 1, 1, 0, 1, 1]
    assert prefix(len(displayed)) == displayed

    # Check the defining digit extraction and the 2-3-block recurrence on a
    # substantial exact finite range, including every residue class tested.
    limit = 100_000
    for n in range(limit):
        digits = sesquinary_digits(n)
        assert all(0 <= digit <= 2 for digit in digits)
        assert digits[-1] != 0 or n == 0
        assert sesquinary_value(digits) == n
        assert x(3 * n) == x(2 * n)
        assert x(3 * n + 1) == 1 - x(2 * n)
        assert x(3 * n + 2) == x(2 * n)

    sample = prefix(1_000_000)
    counts = {factor: count_factor(sample, factor)
              for factor in ((0, 0), (0, 1), (1, 0), (1, 1))}
    # Exact rational error bounds around the two conjectured values.  These
    # are diagnostic finite checks, not evidence of a limit by themselves.
    total_starts = len(sample) - 1
    assert abs(Fraction(counts[(0, 0)], total_starts) - Fraction(1, 10)) < Fraction(1, 100)
    assert abs(Fraction(counts[(0, 1)], total_starts) - Fraction(4, 10)) < Fraction(1, 100)
    assert abs(Fraction(counts[(1, 1)], total_starts) - Fraction(1, 10)) < Fraction(1, 100)
    assert abs(Fraction(counts[(1, 0)], total_starts) - Fraction(4, 10)) < Fraction(1, 100)

    print("SQ(7)=", "".join(map(str, reversed(sesquinary_digits(7)))))
    print("prefix_25=", "".join(map(str, displayed)))
    print("factor_counts_N=", len(sample), counts)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
