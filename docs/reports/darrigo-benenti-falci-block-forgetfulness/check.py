#!/usr/bin/env python3
"""Exact checks for the D'Arrigo--Benenti--Falci forgetfulness conjecture.

The checks use SymPy expressions and rational arithmetic only.  The choice
tau = log(3), tau_p = tau_c = 1 gives r = exp(-tau/tau_c) = 1/3 exactly;
tau > tau_p is checked symbolically, so the covariance windows do not
overlap in the geometric-tail calculation.
"""

from __future__ import annotations

from fractions import Fraction
import sympy as sp


def self_covariance(tau_p: sp.Expr, tau_c: sp.Expr) -> sp.Expr:
    """Integral of exp(-|v-u|/tau_c) over one transit window."""
    return 2 * (tau_c * tau_p - tau_c**2 * (1 - sp.exp(-tau_p / tau_c)))


def separated_covariance(distance: sp.Expr, tau_p: sp.Expr, tau_c: sp.Expr) -> sp.Expr:
    """Same double integral when distance >= tau_p (the case used below)."""
    return (
        tau_c**2
        * sp.exp(-distance / tau_c)
        * (sp.exp(tau_p / tau_c) - 1)
        * (1 - sp.exp(-tau_p / tau_c))
    )


def boundary_sum(n: int, m: int, idle: int, r: sp.Rational) -> sp.Expr:
    """Cross-block covariance sum after factoring B=(1-r)^2/r."""
    total = sp.Integer(0)
    for b0 in range(m):
        for b1 in range(b0 + 1, m):
            for u in range(n):
                for v in range(n):
                    # Use indices b*(n+L)+u.  For tau=tau_c=1 and
                    # tau_p=1, K_ij = (1-r)^2*r^(distance-1).
                    distance = (b1 - b0) * (n + idle) + v - u
                    total += (1 - r) ** 2 * r ** (distance - 1)
    return sp.factor(total)


def check_fidelity() -> None:
    # The source's one-use factor is g = exp(-2 K_00), because a sigma_z
    # eigenvalue difference is 2.  For |+><+|, F(|+>, Phi(|+><+|))=(1+g)/2.
    tau_p = sp.Integer(1)
    tau_c = sp.Integer(1)
    k00 = self_covariance(tau_p, tau_c)
    g = sp.exp(-2 * k00)
    fidelity = sp.simplify((1 + g) / 2)
    assert sp.simplify(k00 - 2 * sp.exp(-1)) == 0
    assert sp.simplify(fidelity - (1 + sp.exp(-4 * sp.exp(-1))) / 2) == 0
    print("fidelity:", fidelity)


def check_covariance_cut() -> None:
    # Exact Lorentzian parameters: tau=tau_c*log(3), tau_p=tau_c=1.
    # Thus r=exp(-tau/tau_c)=1/3, B=(1-r)^2/r=4/3, and log(3)>1.
    tau = sp.log(3)
    assert sp.ask(sp.Q.positive(tau - 1)) is True
    r = sp.Rational(1, 3)
    B = (1 - r) ** 2 / r
    assert B == sp.Rational(4, 3)

    checked = 0
    for n in range(1, 5):
        for m in range(2, 6):
            for idle in range(0, 5):
                cross = boundary_sum(n, m, idle, r)
                bound = sp.Integer(m - 1) * r**idle
                # This is exactly Lemma 4 after substituting B and
                # (1-r)^(-2); all quantities remain rational.
                assert sp.simplify(bound - cross) >= 0
                checked += 1
    print("boundary_cases:", checked)


def check_source_constants() -> None:
    # Eq. (38) has h=4 lambda^2 g^2 tau_c^2.  Keep it symbolic and verify
    # the displayed Lorentzian factor is bounded by h*r^L for L=1,2,3.
    lam = sp.Integer(1)
    tau_p = sp.Integer(1)
    tau_c = sp.Integer(1)
    g = sp.exp(-4 * sp.exp(-1))
    h = 4 * lam**2 * g**2 * tau_c**2
    r = sp.Rational(1, 3)
    transit = (1 - sp.exp(-tau_p / tau_c)) ** 2
    for idle in (1, 2, 3):
        rhs = h * r**idle
        displayed = h * transit * r**idle
        assert sp.simplify((rhs - displayed) - h * r**idle * (1 - transit)) == 0
        assert sp.ask(sp.Q.positive(rhs - displayed)) is True


if __name__ == "__main__":
    check_fidelity()
    check_covariance_cut()
    check_source_constants()
    print("ALL_CHECKS_OK")
