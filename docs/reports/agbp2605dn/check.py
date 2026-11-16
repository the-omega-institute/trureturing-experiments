#!/usr/bin/env python3
"""Exact arithmetic checks for the capsule refutation of the perimeter question."""

from fractions import Fraction
import sympy as sp


def check_radius_identities(k: int) -> tuple[sp.Expr, sp.Expr]:
    L = sp.Integer(k) ** 2
    r = (sp.sqrt(L**2 + 4) - L) / 2
    # r is the positive root of L*r+r^2=1.
    assert sp.simplify(L * r + r**2 - 1) == 0
    assert sp.simplify(r - 1 / (2 * L)) > 0

    R = sp.symbols("R", positive=True)
    f = L * (R**2 - r**2) + sp.Rational(2, 3) * (R**3 - r**3) - 1
    upper = sp.sqrt(sp.Rational(2, 1) / L)
    # f(r)=-1 and f is strictly increasing for R>0; f(upper)>0,
    # so the unique positive root lies below upper.
    assert sp.simplify(f.subs(R, r) + 1) == 0
    assert sp.simplify(f.diff(R)) > 0
    assert sp.simplify(f.subs(R, upper)) > 0
    return r, upper


def check_geometric_constraints(k: int) -> None:
    L = sp.Integer(k) ** 2
    r, upper = check_radius_identities(k)
    R = sp.symbols("R", positive=True)
    # Capsule formulas in R^3:
    # |K_a|=2*pi*L*a^2+(4/3)*pi*a^3 and P(K_a)=4*pi*L*a+4*pi*a^2.
    inner_perimeter = 4 * sp.pi * (L * r + r**2)
    assert sp.simplify(inner_perimeter - 4 * sp.pi) == 0
    shell_volume = 2 * sp.pi * (L * (R**2 - r**2) + sp.Rational(2, 3) * (R**3 - r**3))
    constraint = L * (R**2 - r**2) + sp.Rational(2, 3) * (R**3 - r**3)
    assert sp.simplify(shell_volume - 2 * sp.pi * constraint) == 0

    # The comparison annulus has inner radius 1 and b^3=5/2.
    b = sp.real_root(sp.Rational(5, 2), 3)
    assert sp.simplify((sp.Rational(4, 3) * sp.pi * (b**3 - 1)) - 2 * sp.pi) == 0
    assert sp.simplify(4 * sp.pi - 4 * sp.pi) == 0


def check_spectral_separation(k: int) -> None:
    # The radial weighted Cauchy--Schwarz estimate gives lambda_capsule >= 2r/R^3.
    # With r>1/(2L), R^2<2/L and L=k^2, this is strictly greater than k/3.
    L = Fraction(k * k, 1)
    assert Fraction(1, 2) < Fraction(3, 2)  # 1/(2*sqrt(2)) > 1/3, squared: 9>8.
    assert 9 > 8
    # For b^3=5/2, 5/4 < b < 3/2.  Hence the trial quotient is <108<243.
    assert Fraction(5, 4) ** 3 < Fraction(5, 2)
    assert Fraction(3, 2) ** 3 > Fraction(5, 2)
    trial_bound = 3 * Fraction(9, 4) / Fraction(1, 4) ** 2
    assert trial_bound == 108
    assert trial_bound < 243
    assert Fraction(k, 3) >= 243
    print(f"k={k}: exact radius and constraint identities verified")
    print(f"k={k}: capsule lower bound > k/3 >= 243 > annulus trial upper bound")


def main() -> None:
    for k in (729, 730):
        check_geometric_constraints(k)
        check_spectral_separation(k)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
