#!/usr/bin/env python3
"""Exact checks for the Tarabunga et al. natural-orbital P.4 refutation.

All evidence is rational/integer.  Logarithms occur only in explanatory output;
the decisive Renyi inequality is compared after exponentiation as integers.
"""

from fractions import Fraction
from math import gcd


def h3_argument(p: Fraction) -> Fraction:
    """The positive argument of h_3(p) = -1/2 log_2(argument)."""
    return 1 - 3 * p + 3 * p * p


def assert_equal(label, left, right):
    if left != right:
        raise AssertionError(f"{label}: {left!r} != {right!r}")
    print(f"{label}: {left}")


def main():
    # Source fidelity: for m >= 4 a quadratic Majorana monomial changes at most
    # two occupations, so it has zero matrix element between |0^m> and |1^m>.
    # Hence G_m(p,phi) has Gamma=(1-2p)J_m and every normal-form basis gives
    # the two probabilities 1-p,p.
    m = 6
    p = Fraction(1, 20)
    assert m >= 4 and 0 < p < Fraction(1, 2)
    gamma_scalar = 1 - 2 * p
    assert_equal("G_6 covariance scalar", gamma_scalar, Fraction(9, 10))
    assert_equal("G_6 binary probability sum", (1 - p) + p, Fraction(1))
    assert_equal("G_6 h3 argument", h3_argument(p), Fraction(343, 400))

    # Exact six-mode Gaussian circuit witness from the preregistration.
    # The four squared amplitudes are 57,19,1,3 divided by 80.
    amplitudes_sq = [Fraction(57, 80), Fraction(19, 80), Fraction(1, 80), Fraction(3, 80)]
    assert_equal("witness norm", sum(amplitudes_sq), Fraction(1))
    assert all(x.denominator == 80 for x in amplitudes_sq)
    assert_equal("outcome q0", amplitudes_sq[0] + amplitudes_sq[2], Fraction(29, 40))
    assert_equal("outcome q1", amplitudes_sq[1] + amplitudes_sq[3], Fraction(11, 40))

    q0, q1 = Fraction(29, 40), Fraction(11, 40)
    p0, p1 = Fraction(1, 58), Fraction(3, 22)
    assert_equal("branch-parameter average", q0 * p0 + q1 * p1, p)
    assert_equal("h3 branch 0 argument", h3_argument(p0), Fraction(3193, 3364))
    assert_equal("h3 branch 1 argument", h3_argument(p1), Fraction(313, 484))
    assert 0 < p0 < p < p1 < Fraction(1, 6)

    # Strict Jensen gap, written without floating point or an unverified log.
    # Since h_3(x)=-1/2 log_2(1-3x+3x^2), the desired average increase is
    # equivalent to this integer inequality after raising to the 40th power.
    lhs = 3193**29 * 313**11 * 400**40
    rhs = 3364**29 * 484**11 * 343**40
    if not lhs < rhs:
        raise AssertionError("exact Renyi-3 certificate failed")
    print("exact Renyi-3 certificate: lhs < rhs")
    print(f"certificate bit lengths: lhs={lhs.bit_length()}, rhs={rhs.bit_length()}")

    # Family formulas: q0 p0 + q1 p1 = p for a second exact parameter.
    p_family = Fraction(1, 32)
    q0f = (3 - 2 * p_family) / 4
    q1f = (1 + 2 * p_family) / 4
    p0f = p_family / (3 - 2 * p_family)
    p1f = 3 * p_family / (1 + 2 * p_family)
    assert_equal("family weighted parameter", q0f * p0f + q1f * p1f, p_family)
    assert 0 < p0f < p_family < p1f < Fraction(1, 6)
    print("family parameter ordering: 0 < p0 < p < p1 < 1/6")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
