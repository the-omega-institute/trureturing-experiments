#!/usr/bin/env python3
"""Exact checks for Oudra's conjecture in OEIS A248186, revision 18.

Only integer and fractions.Fraction arithmetic is used. The minimum is
computed from the literal finite sum, independently of the proposed formula.
The real cube root is isolated between rationals by integer cube comparisons.
"""

from fractions import Fraction as F


# OEIS A248186 revision 18, DATA (%S, %T, %U); offset 1.
SOURCE_TERMS = (
    1, 1, 1, 1, 2, 3, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 10, 11,
    12, 12, 13, 14, 14, 15, 16, 17, 17, 18, 19, 19, 20, 21, 21, 22,
    23, 23, 24, 25, 26, 26, 27, 28, 28, 29, 30, 30, 31, 32, 32, 33, 34,
    35, 35, 36, 37, 37, 38, 39, 39, 40, 41, 41, 42, 43, 44, 44, 45, 46,
)
LIMIT = 10000


def integer_cube_root(x):
    """Return floor(cube_root(x)), using only integer comparisons."""
    assert x >= 0
    lo, hi = 0, 1 << ((x.bit_length() + 2) // 3)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** 3 <= x:
            lo = mid
        else:
            hi = mid
    assert lo ** 3 <= x < (lo + 1) ** 3
    return lo


def root_interval():
    denominator = 10 ** 80
    numerator = integer_cube_root(denominator ** 3 // 3)
    lo, hi = F(numerator, denominator), F(numerator + 1, denominator)
    assert lo ** 3 < F(1, 3) < hi ** 3
    return lo, hi


def conjectured_floor(n, lo, hi):
    assert n >= 3 and 0 < lo < hi
    lower = lo * n + lo ** 2 / n - 1
    upper = hi * n + hi ** 2 / n - 1
    left = lower.numerator // lower.denominator
    right = upper.numerator // upper.denominator
    assert left == right, ("root interval does not certify floor", n)
    return left


def main():
    # Direct sums, with no telescoping formula in their construction.
    residuals = [F(1, 18)]
    for k in range(1, LIMIT + 1):
        residuals.append(residuals[-1] - F(1, k * (k + 1) * (k + 2) * (k + 3)))
        assert residuals[-1] > 0
        assert residuals[-1] < residuals[-2]
    assert all(residuals[k] == F(1, 3 * (k + 1) * (k + 2) * (k + 3))
               for k in range(LIMIT + 1))
    print(f"TELESCOPING_CHECK: literal sums agree with residual formula for 0 <= k <= {LIMIT}")

    # Source EXAMPLE: a(5)=2, with both comparisons made exactly.
    assert residuals[2] == F(1, 180) < F(1, 125) < F(1, 72) == residuals[1]
    print("FIDELITY_EXAMPLE: a(5)=2; 1/180 < 1/125 < 1/72")

    # The source's five maximum-n examples, without decimal approximations.
    max_n_examples = ((1, 4), (2, 5), (3, 7), (4, 8), (5, 10))
    for k, n in max_n_examples:
        assert residuals[k] < F(1, n ** 3)
        assert residuals[k] >= F(1, (n + 1) ** 3)
    print("FIDELITY_TABLE: (k, largest n) = (1,4), (2,5), (3,7), (4,8), (5,10)")

    lo, hi = root_interval()
    actual = []
    k = 1
    for n in range(1, LIMIT + 1):
        threshold = F(1, n ** 3)
        while residuals[k] >= threshold:
            k += 1
        # Positivity and every smaller positive candidate are checked directly.
        assert k >= 1 and residuals[k] < threshold
        if k > 1:
            assert residuals[k - 1] >= threshold
        actual.append(k)
        if n >= 3:
            assert k == conjectured_floor(n, lo, hi), (n, k)
    assert tuple(actual[:len(SOURCE_TERMS)]) == SOURCE_TERMS
    print(f"FIDELITY_DATA: all {len(SOURCE_TERMS)} displayed source terms reproduced")
    print(f"CONJECTURE_CHECK: exact agreement for every 3 <= n <= {LIMIT} ({LIMIT - 2} cases)")
    print("ROOT_CERTIFICATE: rational interval of width 1/10^80; both endpoint floors agree in every case")
    print(f"ENDPOINTS: a(1)={actual[0]}, a(2)={actual[1]}, a(3)={actual[2]}, a({LIMIT})={actual[-1]}")
    print("SCOPE: finite exact checks support the outline; no unbounded or Lean proof is claimed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
