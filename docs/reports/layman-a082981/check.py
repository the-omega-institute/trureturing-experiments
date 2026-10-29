#!/usr/bin/env python3
"""Exact finite checks for OEIS A082981 (John W. Layman, 2003)."""
from fractions import Fraction

OEIS_PREFIX = [
    1, 2, 3, 4, 9, 14, 19, 24, 53, 82, 111, 140, 309, 478, 647, 816,
    1801, 2786, 3771, 4756, 10497, 16238, 21979, 27720, 61181, 94642,
    128103, 161564, 356589, 551614, 746639, 941664, 2078353, 3215042,
    4351731, 5488420, 12113529, 18738638, 25363747, 31988856, 70602821,
]


def active(x: int, y: int) -> bool:
    # The source condition is 4*y <= 5*x, including equality.
    return 5 * x >= 4 * y


def insert_stage(stage: list[int]) -> list[int]:
    out = [stage[0]]
    for x, y in zip(stage, stage[1:]):
        if active(x, y):
            out.append(x + y)
        out.append(y)
    return out


def finite_process(rounds: int = 250, prefix_capacity: int = 3000) -> list[int]:
    """Run the simultaneous rule, retaining a large exact prefix.

    A retained prefix through 3000 entries is enough for the requested 200
    terms in this finite reproduction. The check also compares the last two
    rounds, so the reported prefix has stopped changing within the run.
    """
    stage = [1, 1]
    for _ in range(rounds):
        stage = insert_stage(stage)[:prefix_capacity]
    return stage


def recurrence_terms(count: int) -> list[int]:
    terms = [1, 2]
    for i in range(count - 2):
        ratio = Fraction(terms[-2], terms[-1])
        q = 1 + (Fraction(5, 4) + ratio).numerator // (Fraction(5, 4) + ratio).denominator
        terms.append(q * terms[-1] - terms[-2])
    return terms


def main() -> None:
    stage = [1, 1]
    previous = None
    for _ in range(250):
        previous = stage
        stage = insert_stage(stage)[:3000]
    assert stage[:200] == previous[:200], "the first 200 terms changed in the final round"

    generated = stage[:200]
    assert generated[: len(OEIS_PREFIX)] == OEIS_PREFIX

    terms = recurrence_terms(200)
    assert generated == terms
    assert terms[:8] == [1, 2, 3, 4, 9, 14, 19, 24]
    assert all(terms[n + 8] == 6 * terms[n + 4] - terms[n] for n in range(192))

    # The four exact ratio phases are I0, I1, I2, I3, with q=(2,2,3,2).
    phases = []
    for i in range(12):
        r = Fraction(terms[i], terms[i + 1])
        phases.append((r, 1 + (Fraction(5, 4) + r).numerator // (Fraction(5, 4) + r).denominator))
    assert [q for _, q in phases] == [2, 2, 3, 2] * 3

    # Exact return-map arithmetic from the proposed proof.
    B = ((-5, 7), (-8, 11))
    B2 = tuple(tuple(sum(B[i][k] * B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
    assert B2 == ((-31, 42), (-48, 65))
    assert tuple(tuple(B2[i][j] - 6 * B[i][j] + (1 if i == j else 0) for j in range(2)) for i in range(2)) == ((0, 0), (0, 0))

    print("rounds=250 retained_prefix=3000 stable_terms=200")
    print("oeis_revision_22_prefix_match=41")
    print("lag_four_identities_checked=192")
    print("phase_coefficients=2,2,3,2 repeated")
    print("B2_minus_6B_plus_I=zero_matrix")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
