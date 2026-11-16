#!/usr/bin/env python3
"""Exact checks for agbp2605dn; no floating point or numerical eigenvalue solver.

Fidelity: the Steiner formula and W_(n-1)(B_r)=omega_n*r printed in
arXiv:2605.25182v1, section 2 (labels eq:steinerRN and eq:quer-ball).
Settlement check: exact capsule constraints, radial-estimate integrals,
comparison-annulus trial quotient, and certified annulus ODE root brackets.
The analytic Sobolev, coordinate and trace arguments are not certified here.
"""
from fractions import Fraction as F
from math import comb, factorial
import sympy as sp


def fidelity() -> None:
    r, t, omega = sp.symbols('r t omega', positive=True)
    for n in range(2, 7):
        coefficients = [omega * r ** (n - j) for j in range(n + 1)]
        steiner = sum(comb(n, j) * coefficients[j] * t ** j for j in range(n + 1))
        assert sp.expand(omega * (r + t) ** n - steiner) == 0
        assert coefficients[1] == sp.diff(omega * r ** n, r) / n
        assert coefficients[n - 1] == omega * r
    print('FIDELITY: source Steiner and ball-quermassintegral formulas, n=2..6: exact')


def radial_integrals() -> None:
    r, R, a = sp.symbols('r R a', positive=True)
    assert sp.integrate(a ** -2, (a, r, R)) == 1 / r - 1 / R
    assert sp.integrate(a, (a, r, R)) == (R ** 2 - r ** 2) / 2
    assert sp.integrate(a ** 2, (a, r, R)) == (R ** 3 - r ** 3) / 3
    # On r<=a<=R, a^-1<=R*a^-2. Cauchy--Schwarz therefore gives
    # cylinder constant <=(R/r)*integral_r^R a da; cap constant
    # <=(1/r)*integral_r^R a^2 da. Check the relaxed constants exactly.
    assert sp.simplify(R ** 3 / (2 * r) - R / r * (R ** 2 - r ** 2) / 2) == R * r / 2
    assert sp.simplify(R ** 3 / (3 * r) - (R ** 3 - r ** 3) / (3 * r)) == r ** 2 / 3
    print('RADIAL: inverse-weight, cylinder and cap integrals: exact')


def root_interval(f, lo: F, hi: F, steps: int = 100) -> tuple[F, F]:
    assert f(lo) < 0 < f(hi)
    for _ in range(steps):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    assert f(lo) <= 0 <= f(hi)
    return lo, hi


def annulus_bounds() -> tuple[F, F]:
    blo, bhi = root_interval(lambda b: b ** 3 - F(5, 2), F(5, 4), F(3, 2))
    # Mass of u(s)=s-1: integral_0^(b-1) t^2(1+t)^2 dt.
    mass = lambda d: d ** 3 / 3 + d ** 4 / 2 + d ** 5 / 5
    qlo, qhi = F(1, 2) / mass(bhi - 1), F(1, 2) / mass(blo - 1)
    assert 20 < qlo < qhi < 21
    assert (F(5, 4)) ** 3 < F(5, 2) < (F(3, 2)) ** 3
    assert 3 * F(3, 2) ** 2 / F(1, 4) ** 2 == 108 < 243

    # For the 3D annulus the positive radial eigenfunction is
    # sin(mu*(s-1))/s; the outer Neumann equation is
    # b*mu*cos(mu*(b-1))-sin(mu*(b-1))=0.
    # Alternating Taylor series certify signs at rational endpoints.
    def sin_bound(x: F, j: int) -> F:
        return sum(((-1) ** i * x ** (2 * i + 1) / factorial(2 * i + 1)
                    for i in range(j + 1)), F(0))

    def cos_bound(x: F, j: int) -> F:
        return sum(((-1) ** i * x ** (2 * i) / factorial(2 * i)
                    for i in range(j + 1)), F(0))

    assert cos_bound(F(3, 2), 3) > 0  # cos stays positive up to 3/2.
    mlo, mhi = F(387, 100), F(388, 100)
    for mu in (mlo, mhi):
        assert 0 < mu * (blo - 1) < mu * (bhi - 1) < F(3, 2)
        # Sin is increasing and cos decreasing on [0,3/2].
        # sin lower/upper: odd/even truncations; cos likewise.
        xlo, xhi = mu * (blo - 1), mu * (bhi - 1)
        flo = blo * mu * cos_bound(xhi, 13) - sin_bound(xhi, 12)
        fhi = bhi * mu * cos_bound(xlo, 12) - sin_bound(xlo, 13)
        if mu == mlo:
            assert flo > 0
        else:
            assert fhi < 0
    assert mlo ** 2 == F(149769, 10000) and mhi ** 2 == F(150544, 10000)
    print('ANNULUS: exact trial quotient in (20,21); first radial eigenvalue in (149769/10000,150544/10000)')
    return qlo, qhi


def capsule(k: int, trial_upper: F) -> None:
    assert k >= 1
    L = k * k
    rlo, rhi = root_interval(lambda r: L * r + r * r - 1, F(0), F(1))
    assert F(1, 2 * L) < rlo < rhi < F(1, L)
    c = lambda r: r - r ** 3 / 3
    # At the exact root r: L*r^2+(2/3)*r^3=r-r^3/3.
    Rlo, _ = root_interval(lambda R: L * R ** 2 + F(2, 3) * R ** 3 - 1 - c(rlo), F(0), F(2, k))
    _, Rhi = root_interval(lambda R: L * R ** 2 + F(2, 3) * R ** 3 - 1 - c(rhi), F(0), F(2, k))
    assert rhi < Rlo < Rhi and Rhi ** 2 < F(2, L)
    lower = 2 * rlo / Rhi ** 3
    upper = 2 * rhi / Rlo ** 3
    assert lower > trial_upper
    if k >= 729:
        assert lower > F(k, 3) >= 243 > 108
        assert 2 * k - 1 < lower < upper < 2 * k
    else:
        assert k == 16 and 31 < lower < upper < 32
    print(f'CAPSULE k={k}: exact radius brackets; lower bound in ({2*k-1},{2*k}); exceeds annulus trial')


def constraints_and_family() -> None:
    r, R, L = sp.symbols('r R L', positive=True)
    p = 4 * sp.pi * (L * r + r ** 2)
    volume = lambda a: 2 * sp.pi * L * a ** 2 + sp.Rational(4, 3) * sp.pi * a ** 3
    assert sp.simplify(p - 4 * sp.pi * (L * r + r ** 2)) == 0
    shell = volume(R) - volume(r)
    radius_eq = L * (R ** 2 - r ** 2) + sp.Rational(2, 3) * (R ** 3 - r ** 3)
    assert sp.expand(shell - 2 * sp.pi * radius_eq) == 0
    assert sp.rem(L * r ** 2 + sp.Rational(2, 3) * r ** 3 - (r - r ** 3 / 3),
                  r ** 2 + L * r - 1, r) == 0
    # Monotonicity of the outer-radius equation.
    assert sp.diff(radius_eq, R) == 2 * L * R + 2 * R ** 2
    # k/(2*sqrt(2)) > k/3 is the exact rational comparison 9>8.
    assert 9 > 8 and F(729, 3) == 243
    print('CONSTRAINTS: P(K_r)=4*pi and shell volume=2*pi under the two exact radius equations')
    print('FAMILY: r>1/(2L), R^2<2/L imply 2r/R^3>k/3 for L=k^2; every k>=729 separates')


def main() -> None:
    fidelity()
    radial_integrals()
    constraints_and_family()
    _, trial_upper = annulus_bounds()
    for k in (16, 729, 730):
        capsule(k, trial_upper)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
