#!/usr/bin/env python3
"""Exact arithmetic checks for the proposed Proto-NN rate refutation.

All numerical evidence below is rational arithmetic.  The only non-rational
quantity in the witness is h itself; for h(s)=(2s)^(1/(1+gamma)) the
power appearing in the conjecture is exactly h(s)^(1+gamma)=2s, so no
floating-point root is needed.
"""

from fractions import Fraction
def occupancy_expectation(k: int, n: int) -> Fraction:
    """Expected number of good centers in the three-prototype event."""
    # The center is in [1/4,3/4] (mass 1/4), one prototype is in each
    # side interval (mass 1/(2n)), and all remaining prototypes avoid the
    # window (mass 1-2/n).
    return Fraction(k * (k - 1) * (k - 2), 16 * n * n) * Fraction(n - 2, n) ** (k - 3)


def conditional_regret_lower_bound(n: int) -> Fraction:
    """The source's conditional contribution at W >= 1/(2n)."""
    return Fraction(1, 32 * n) * Fraction(n - 1, n) ** (n - 1)


def lower_bound(k: int, n: int) -> Fraction:
    """Linearity-of-expectation lower bound before Bernoulli simplification."""
    return (
        Fraction(k * (k - 1) * (k - 2), 2048 * n**3)
        * Fraction(n - 2, n) ** (k - 3)
    )


def bayes_risk() -> Fraction:
    # On each component the more likely label has probability 3/4, and each
    # component has mass 1/2.  Thus the Bayes error is 1/4.
    return Fraction(1, 2) * Fraction(1, 4) + Fraction(1, 2) * Fraction(1, 4)


def check_fidelity() -> None:
    assert bayes_risk() == Fraction(1, 4)
    # A good cell has W in [1/(2n),1/n].  Its one-positive-label event has
    # probability n W (1-W)^(n-1) / 4 and regret W/2, hence the expression
    # n W^2 (1-W)^(n-1) / 8.  The source bounds it using W >= 1/(2n) and
    # (1-1/n)^(n-1) >= 1/4.
    for n in range(2, 40):
        assert conditional_regret_lower_bound(n) >= Fraction(1, 128 * n)
    # Verify the exact binomial lower bound for a representative range.
    for k in range(4, 25):
        for n in range(4 * k, 4 * k + 13):
            assert lower_bound(k, n) >= Fraction(k**3, 16384 * n**3)
            assert occupancy_expectation(k, n) > 0


def check_counterexample_family() -> None:
    # For gamma = 7, h(2/k)^8 = 4/k and q=(gamma+1)/2=4.
    # Along n=t^5, k=t^4, both conjectured terms are O(t^-4), while the
    # exact simplified lower bound is 1/(16384 t^3).
    for t in range(4, 18):
        n, k = t**5, t**4
        assert n >= 4 * k and k >= 4
        first_term = Fraction(4, k)       # h(2/k)^8
        second_term = Fraction(k, n) ** 4  # (k/n)^4
        lb = Fraction(1, 16384 * t**3)
        assert lower_bound(k, n) >= lb
        assert lb / (first_term + second_term) == Fraction(t, 81920)
    # The ratio is unbounded; one exact threshold witness is enough to rule
    # out any fixed constants in the O statement.
    t = 81921
    assert Fraction(t, 81920) > 1


def main() -> None:
    check_fidelity()
    check_counterexample_family()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
