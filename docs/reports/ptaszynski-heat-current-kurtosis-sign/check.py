#!/usr/bin/env python3
"""Exact checks for the thermal heat-current kurtosis candidate.

All rates, determinants, derivatives, and sign decisions stay in SymPy's exact
integer/rational domain.  The counting field is represented by x = exp(lambda)
until the final implicit-series calculation.
"""

from fractions import Fraction

import sympy as sp


def tilted_matrix(energies, transitions, x):
    """Build W^h with a transition j -> i recorded by heat q to the bath."""
    n = len(energies)
    w = sp.zeros(n)
    for i, j, rate, q in transitions:
        w[i, j] += sp.sympify(rate) * x**q
        w[j, j] -= sp.sympify(rate)
    return w


def implicit_cumulants(det, z, x):
    """Solve det(chi(lambda) I - W(lambda)) = 0 through order four."""
    lam = sp.symbols("lambda")
    cs = sp.symbols("c1:5")
    chi = sum(cs[n - 1] * lam**n / sp.factorial(n) for n in range(1, 5))
    series = sp.series(det.subs(x, sp.exp(lam)).subs(z, chi), lam, 0, 5)
    equations = sp.expand(series.removeO())
    solved = {}
    for n, c in enumerate(cs, 1):
        coefficient = sp.expand(equations.subs(solved)).coeff(lam, n)
        solved[c] = sp.factor(sp.solve(coefficient, c)[0])
    return tuple(solved[c] for c in cs)


def family_transitions(R):
    # (L0,H0,L1,H1,L2,H2), with the hot edges followed by cold edges.
    transitions = []
    for L, H, k in ((0, 1, 1), (2, 3, 6), (4, 5, 6)):
        transitions += [(H, L, k, 0), (L, H, 2 * k, 0)]
    for H, L in ((1, 2), (3, 4), (5, 0)):
        transitions += [(L, H, R, 1), (H, L, 1, -1)]
    return transitions


def family_check():
    R, z, x = sp.symbols("R z x", positive=True)
    energies = (0, 1, 0, 1, 0, 1)
    transitions = family_transitions(R)

    # Local detailed balance: hot upward/downward = 1/2, cold = 1/R.
    assert sp.Rational(1, 2) == sp.exp(-sp.log(2))
    assert sp.simplify(R**-1 - sp.exp(-sp.log(R))) == 0
    assert all(rate > 0 for _, _, rate, _ in transitions if sp.sympify(rate).is_number)

    W = tilted_matrix(energies, transitions, x)
    det = sp.factor((z * sp.eye(6) - W).det())
    q = R**2 + 2 * R + 4
    coeff = {
        1: 1488 + 936 * R + 468 * R**2 + 48 * R**3,
        2: 2738 + 1997 * R + 388 * R**2 + 13 * R**3,
        3: 2072 + 879 * R + 68 * R**2 + R**3,
        4: 539 + 97 * R + 3 * R**2,
        5: 42 + 3 * R,
        6: 1,
    }
    P = sum(coeff[j] * z**j for j in range(1, 7))
    expected_det = P - 36 * R**3 * (x**3 - 1) - 288 * (x**-3 - 1)
    assert sp.factor(det - expected_det) == 0

    F, B = 36 * R**3, sp.Integer(288)
    d1, d2, d3, d4 = (coeff[j] for j in (1, 2, 3, 4))
    c1 = 3 * (F - B) / d1
    c2 = (9 * (F + B) - 2 * d2 * c1**2) / d1
    c3 = (27 * (F - B) - 6 * d2 * c1 * c2 - 6 * d3 * c1**3) / d1
    c4 = (
        81 * (F + B)
        - d2 * (8 * c1 * c3 + 6 * c2**2)
        - 36 * d3 * c1**2 * c2
        - 24 * coeff[4] * c1**4
    ) / d1
    r8 = {R: 8}
    assert (sp.factor(c1.subs(r8)), sp.factor(c2.subs(r8))) == (
        sp.Rational(6, 7),
        sp.Rational(3581, 2401),
    )
    assert sp.factor(c3.subs(r8)) == sp.Rational(1349819, 1647086)
    assert sp.factor(c4.subs(r8)) == sp.Rational(-113056891, 1129900996)
    assert sp.factor((c4 / c2).subs(r8)) == sp.Rational(-113056891, 1685204276)

    # The source's half-line certificate: c4 has a negative polynomial in R-8.
    a = (
        7690632602871168,
        19461033025247808,
        15323367935019744,
        6343945602823440,
        1649604944057256,
        293504823956556,
        37474457340882,
        3526498268247,
        247179395436,
        12798157749,
        471833962,
        11351565,
        145432,
        491,
    )
    y = sp.symbols("y")
    p = -sum(sp.Integer(v) * y**j for j, v in enumerate(a))
    denominator = 8 * (4 * R + 31) ** 7 * q**3
    assert sp.factor(c4 - 243 * p.subs(y, R - 8) / denominator) == 0
    assert all(v > 0 for v in a)
    assert sp.factor(c2.subs(R, 8)) > 0
    print("family R=8: c2 =", sp.factor(c2.subs(R, 8)))
    print("family R=8: c4/c2 =", sp.factor((c4 / c2).subs(R, 8)))
    print("family R>=8: c4 polynomial certificate has", len(a), "positive coefficients")


def second_witness_check():
    z, x = sp.symbols("z x")
    energies = (0, 5, 1, 6, 2, 7)
    transitions = []
    for L, H, forward, reverse in ((0, 1, 4, 8), (2, 3, 24, 48), (4, 5, 24, 48)):
        transitions += [(H, L, forward, 0), (L, H, reverse, 0)]
    for H, L, forward, reverse in ((1, 2, 128, 8), (3, 4, 128, 8), (5, 0, 128, 1)):
        transitions += [(L, H, forward, energies[H] - energies[L])]
        transitions += [(H, L, reverse, energies[L] - energies[H])]
    # Hot ratios are exp(-(1/5) log 2 * gap), cold ratios are exp(-log 2 * gap).
    assert all(sp.Rational(r, f) == sp.Rational(1, 2) for _, _, f, r in ((0, 0, 8, 4), (0, 0, 48, 24), (0, 0, 48, 24)))
    W = tilted_matrix(energies, transitions, x)
    det = (z * sp.eye(6) - W).det()
    c1, c2, c3, c4 = implicit_cumulants(det, z, x)
    assert c2 == sp.Rational(9010078778600, 35431447851)
    assert c4 == sp.Rational(-70265274785484432257383595000, 53520935151847443648463233)
    print("six distinct energies: c2 =", c2)
    print("six distinct energies: c4 =", c4)


def fidelity_two_state():
    # A two-state, two-temperature example satisfies the source's Eq. (19).
    z, x = sp.symbols("z x")
    transitions = ((1, 0, 1, 0), (0, 1, 2, 0), (1, 0, 1, 1), (0, 1, 4, -1))
    W = tilted_matrix((0, 1), transitions, x)
    cumulants = implicit_cumulants((z * sp.eye(2) - W).det(), z, x)
    K = sp.factor(cumulants[3] / cumulants[1])
    assert K == sp.Rational(5201, 12032) and K > 0
    print("source Eq. (19) fidelity case: K =", K)


if __name__ == "__main__":
    fidelity_two_state()
    family_check()
    second_witness_check()
    print("ALL_CHECKS_OK")
