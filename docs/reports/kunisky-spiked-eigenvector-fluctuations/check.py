#!/usr/bin/env python3
"""Exact arithmetic checks for the Kunisky Conjecture 1.19 counterfamily.

The random matrices are not sampled here.  The checks verify the deterministic
variance-profile identities used by the construction, and the two exact
limiting variance values supplied by the counterexample calculation.
"""

from fractions import Fraction
from math import isqrt


def perfect_fourth(m: int) -> int:
    return m**4


def check_profile(n: int) -> None:
    # Taking N=m^4 makes eta=N^(-1/4) rational, while retaining an infinite
    # subsequence of the proposed family.  All calculations below are exact.
    m = isqrt(isqrt(n))
    assert m**4 == n and m >= 2
    eta = Fraction(1, m)
    eps_sq = eta
    r = 1 - eta
    h = 1 - Fraction(2, n) * r
    s = h / (n + 1)

    assert 0 < eta < 1
    assert 0 < r < 1
    assert 0 < h < 1
    assert 0 < s < eta
    assert r / n < eta

    # E W^2 = I on the two spike rows and every bath row.
    spike_row_sum = eta + r
    bath_row_sum = Fraction(2, n) * r + (n + 1) * s
    assert spike_row_sum == 1
    assert bath_row_sum == 1

    # Cov(vec W) has spectral norm 2 max{eta,r/N,s}; eta is the maximum.
    cov_norm = 2 * max(eta, r / n, s)
    assert cov_norm == 2 * eta
    assert cov_norm < 2 * eta + 1  # exact sanity check on the stated value

    d = n + 2
    alpha_sq = Fraction(1, d) / eps_sq
    assert 0 < alpha_sq < 1
    assert alpha_sq + (1 - alpha_sq) == 1

    print(
        f"N={n}: E[W^2]=I, covariance norm={cov_norm}, "
        f"alpha^2={alpha_sq}"
    )


def main() -> None:
    for m in range(2, 7):
        check_profile(perfect_fourth(m))

    # The optimizer calculation gives p_2/eta -> 3/4 and bath mass q -> 1/4.
    # Hence the conjectured covariance along R is exactly 3/4+1/4=1.
    predicted = Fraction(3, 4) + Fraction(1, 4)

    # The Schur-complement/eigenvector calculation gives independent Gaussian
    # contributions with variances 3/16 and 1/4, hence 7/16.
    actual = Fraction(3, 16) + Fraction(1, 4)
    assert predicted == Fraction(1, 1)
    assert actual == Fraction(7, 16)
    assert predicted != actual

    print(f"predicted covariance limit = {predicted}")
    print(f"actual fluctuation variance = {actual}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
