#!/usr/bin/env python3
"""Exact small checks for Andersson, arXiv:2609.12859v1.

The checks use the finite binomial moment identity for |z_1+z_2|^(2k)
and an exact order/rationality witness for the proposed F_sigma answer.
No floating point values are used as evidence.
"""

from fractions import Fraction
from math import comb

from sympy import Integer, Rational, sqrt


def even_moment(k: int):
    """E |z_1+z_2|^(2k) on the two-torus, exactly."""
    # Only terms with equal powers of z_1 and conjugate z_1 survive.
    return Integer(comb(2 * k, k))


def phi_norm_sq_at_even_exponent(k: int):
    """The square of ||phi||_{2k}^k, represented exactly."""
    # ||phi||_(2k)^(2k) = C(2k,k), hence ||phi||_(2k)^k = sqrt(C(2k,k)).
    return sqrt(even_moment(k))


def source_norm_identity_check() -> None:
    # Source Lemma 1.2 and the displayed T_{n,m} formula at n=3,m=2,p=2.
    # ||phi||_2 = sqrt(2), ||phi||_4^2 = sqrt(6), so ||T_{3,2}|| = sqrt(3)/2.
    phi2 = sqrt(even_moment(1))
    phi4_squared = phi_norm_sq_at_even_exponent(2)
    t_squared = (phi4_squared / phi2**3) ** 2
    assert t_squared == Rational(3, 4)
    assert phi2**2 == Integer(2)
    assert phi4_squared**2 == Integer(6)


def proposed_f_sigma_answer_fails_monotonicity() -> None:
    # The proposed answer explicitly includes B = Q ∩ [1,∞). Any boundedness
    # exponent set of one common operator is upward closed: a bounded extension
    # on H^p restricts to H^q for q >= p.  Here 2 ∈ B but sqrt(5) ∉ B and
    # 2 < sqrt(5), an exact contradiction.
    p = Fraction(2, 1)
    assert p.denominator == 1
    assert sqrt(5) > Integer(2)
    assert sqrt(5).is_rational is False


def main() -> None:
    source_norm_identity_check()
    proposed_f_sigma_answer_fails_monotonicity()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
