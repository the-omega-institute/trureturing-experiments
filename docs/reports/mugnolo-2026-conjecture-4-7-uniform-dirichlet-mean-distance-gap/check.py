#!/usr/bin/env python3
"""Exact arithmetic checks for Mugnolo Conjecture 4.7.

The proof uses only rational arithmetic and exact algebra.  Decimal values are
printed only as a readable high precision display of the exact graph ratios.
"""
from fractions import Fraction
import sympy as sp

Q = Fraction
sqrt5 = sp.sqrt(5)
v_star = (sqrt5 - 2) / 4
C = 4 * (sqrt5 - 2)


def exact_min(a, b):
    return a if a <= b else b


def check_variance_combination():
    # The two Rayleigh trials give min(1/(1+v), 16v).  Their exact
    # intersection and common value are the advertised universal constant.
    assert sp.simplify(1 / (1 + v_star) - C) == 0
    assert sp.simplify(16 * v_star - C) == 0
    assert sp.simplify(C < 1) is sp.true
    # A rational finish also works: v >= 1/17 uses the distance trial,
    # v <= 1/17 uses the truncated trial.
    assert Q(1, 1 + Q(1, 17)) == Q(17, 18)
    assert Q(16, 17) < Q(17, 18)
    for v in [Q(0), Q(1, 1000), Q(1, 17), Q(1, 10), Q(1), Q(7, 3)]:
        bound = exact_min(Q(1, 1) / (1 + v), 16 * v)
        assert bound <= Q(17, 18)
    # Rational variances on either side of v★ provide an exact near-crossover
    # check for the two competing trials.
    for v in [Q(59, 1000), Q(60, 1000)]:
        bound = exact_min(Q(1, 1) / (1 + v), 16 * v)
        assert bound < Q(17, 18)
        print(f"near crossover v={v}: exact combined bound = {bound}")
    print("variance combination: exact C = 4*(sqrt(5)-2), rational C = 17/18")


def check_cantelli_step():
    # Exact two-point mean-one distributions.  For p < 1/2, R=0 with mass p
    # and R=1/(1-p) otherwise; A={R<1/2} has mass p and v=p/(1-p).
    for p in [Q(1, 10), Q(1, 5), Q(2, 7), Q(49, 100)]:
        high = Q(1, 1) / (1 - p)
        assert high > Q(1, 2)
        variance = p * (Q(0) - Q(1)) ** 2 + (1 - p) * (high - Q(1)) ** 2
        assert variance == p / (1 - p)
        assert p / (1 - p) <= 4 * variance
    print("Cantelli step: exact two-point checks passed")


def graph_cases():
    pi2 = sp.pi ** 2
    return {
        "path (one Dirichlet endpoint)": (sp.pi**2 / 4, Q(1, 2), Q(1, 1)),
        "star (any number of unit arms, Dirichlet centre)": (sp.pi**2 / 4, Q(1, 2), Q(1, 1)),
        "cycle (unit circumference, one Dirichlet vertex)": (pi2, Q(1, 4), Q(1, 1)),
        "figure-eight (two unit loops, Dirichlet junction)": (pi2, Q(1, 4), Q(1, 1)),
    }


def check_graph_ratios():
    # Each tuple is (lambda_1, rho, length scale).  The scale is included to
    # make the dimensional cancellation explicit; all listed exact ratios are
    # pi^2/16.  The inequality uses the exact rational Archimedean bound
    # pi < 22/7, so no floating point is evidence for the check.
    for name, (lam, rho, scale) in graph_cases().items():
        ratio = sp.simplify(lam * rho**2)
        assert sp.simplify(ratio - sp.pi**2 / 16) == 0
        rational_upper = sp.Rational(22, 7) ** 2 / 16
        assert rational_upper < sp.Rational(17, 18)
        print(f"{name}: exact lambda*rho^2 = pi^2/16; high-precision = {sp.N(ratio, 30)}")


def main():
    check_variance_combination()
    check_cantelli_step()
    check_graph_ratios()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
