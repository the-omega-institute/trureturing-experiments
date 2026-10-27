#!/usr/bin/env python3
"""Exact checks for Li's Rényi interpolation-convexity counterfamily.

All evidence is exact: integers, Fraction, and SymPy symbolic identities.
"""
from fractions import Fraction
import sympy as sp


def q(rho: Fraction, lam: Fraction) -> Fraction:
    return lam * lam + (1 - lam) * (1 - lam) + 2 * rho * lam * (1 - lam)


def main() -> None:
    rho = Fraction(-1, 2)
    s = Fraction(1, 16)
    t = Fraction(3, 16)
    m = (s + t) / 2

    # Fidelity: projection of the uniform unit-disk density has mass one.
    x = sp.symbols("x", real=True)
    f = sp.Rational(2, 1) / sp.pi * sp.sqrt(1 - x * x)
    mass = sp.simplify(sp.integrate(f, (x, -1, 1)))
    assert mass == 1, mass

    # Fidelity: the quadratic projection scale in the paper's construction.
    r, z = sp.symbols("r z", real=True)
    q_symbolic = z**2 + (1 - z)**2 + 2 * r * z * (1 - z)
    expected = 1 - 2 * (1 - r) * z * (1 - z)
    assert sp.expand(q_symbolic - expected) == 0

    # Exact Jensen defect at the nominated interior points.
    qs, qt, qm = q(rho, s), q(rho, t), q(rho, m)
    assert (qs, qt, qm) == (Fraction(211, 256), Fraction(139, 256), Fraction(172, 256))
    numerator = qm * qm
    denominator = qs * qt
    assert numerator == Fraction(29584, 65536)
    assert denominator == Fraction(29329, 65536)
    assert numerator > denominator
    defect_ratio = numerator / denominator
    assert defect_ratio == Fraction(29584, 29329)

    # The log-profile defect is (1/4) log(defect_ratio), hence is strictly positive.
    defect = sp.log(sp.Rational(defect_ratio.numerator, defect_ratio.denominator)) / 4
    assert defect.is_positive is True

    # Fidelity plus settlement: the curvature numerator is negative near an endpoint.
    # F'' has the sign of rho + 2(1-rho)*lambda*(1-lambda).
    lam = Fraction(1, 32)
    curvature_numerator = rho + 2 * (1 - rho) * lam * (1 - lam)
    assert curvature_numerator < 0
    assert q(rho, lam) > 0

    print("mass_disk_projection = 1")
    print(f"q(1/16) = {qs}")
    print(f"q(3/16) = {qt}")
    print(f"q(1/8) = {qm}")
    print(f"jensen_ratio = {defect_ratio} > 1")
    print("curvature_negative_at_lambda_1_over_32 = True")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
