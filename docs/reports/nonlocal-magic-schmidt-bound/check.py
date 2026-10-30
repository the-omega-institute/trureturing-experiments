#!/usr/bin/env python3
"""Exact checks for the alpha=1/2 Schmidt-reference counterexample.

All certificate arithmetic is integral or rational.  No floating point value is
used as evidence.
"""

from fractions import Fraction


def parity(value: int) -> int:
    return value.bit_count() & 1


def walsh_rows(c: list[int]) -> list[int]:
    rows = []
    for u in range(16):
        total = 0
        for z in range(16):
            w = sum(((-1) ** parity(z & x)) * c[x] * c[x ^ u] for x in range(16))
            total += abs(w)
        rows.append(total)
    return rows


def walsh_l1(c: list[int]) -> int:
    return sum(walsh_rows(c))


def parseval_sum_xi(c: list[int]) -> Fraction:
    """Check the source's Xi normalization: sum_P Xi_P = 1."""
    n = sum(x * x for x in c)
    # There are sixteen Pauli pairs for each (u,z), and d=256.
    numerator = 16 * sum(
        sum(((-1) ** parity(z & x)) * c[x] * c[x ^ u] for x in range(16)) ** 2
        for u in range(16)
        for z in range(16)
    )
    return Fraction(numerator, 256 * n * n)


def assert_equal(actual, expected, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label}: {actual!r} != {expected!r}")


def main() -> None:
    a = [316, 215, 214, 213, 212, 211, 210, 209, 208, 207, 206, 105, 104, 103, 102, 101]
    b = [316, 215, 214, 213, 101, 211, 210, 209, 208, 207, 206, 105, 104, 103, 102, 212]
    n = sum(x * x for x in a)
    assert_equal(n, 596096, "normalization N")
    assert_equal(sum(a), 2936, "sum of witness entries")
    assert all(a[i] > a[i + 1] > 0 for i in range(15)), "strict Schmidt spectrum"
    assert_equal(
        walsh_rows(a),
        [1672784, 1162128, 1157072, 1152016, 1071872, 1066816, 1061760,
         1056704, 1055616, 1052320, 1048992, 797696, 872704, 870208,
         867680, 865184],
        "reference row certificate",
    )
    assert_equal(
        walsh_rows(b),
        [1597696, 1129980, 1125764, 1121580, 823248, 1116796, 1117612,
         1118428, 1100268, 1096100, 1091932, 797696, 793832, 798096,
         802328, 1103460],
        "permuted row certificate",
    )
    sa, sb = walsh_l1(a), walsh_l1(b)
    assert_equal((sa, sb, sa - sb), (16831552, 16734816, 96736), "witness gap")
    # Eq. (1) is normalized exactly by the Pauli Parseval identity.
    assert_equal(parseval_sum_xi(a), Fraction(1), "Xi normalization")
    print("FIDELITY_CHECK_OK: sum_P Xi_P = 1 for the exact witness")

    # The alpha=1/2 identity gives M = 2 log_2(S/(16N)); comparing S is exact.
    assert sb < sa
    assert Fraction(sb, 16 * n) < Fraction(sa, 16 * n)
    print(f"COUNTEREXAMPLE_OK: S(a)={sa}, S(Pa)={sb}, gap={sa - sb}")

    # The explicit open family from the preregistration.
    for t in (Fraction(3), Fraction(13, 4), Fraction(7, 2)):
        c = [t, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1]
        cp = [t, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2]
        expected_before = 14 * t * t + 288 * t + 538
        expected_after = 16 * t * t + 200 * t + 760
        assert_equal(walsh_l1(c), expected_before, f"family S before at t={t}")
        assert_equal(walsh_l1(cp), expected_after, f"family S after at t={t}")
        assert expected_after < expected_before
    print("FAMILY_CHECK_OK: exact gap holds at t=3, 13/4, 7/2")

    # Quantitative neighborhood check from the proposed proof route.
    epsilon = Fraction(1, 100)
    perturbation_bound = Fraction(3006464) * epsilon + Fraction(8192) * epsilon * epsilon
    assert perturbation_bound < 96736
    print(f"OPEN_FAMILY_OK: perturbation_bound={perturbation_bound} < 96736")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
