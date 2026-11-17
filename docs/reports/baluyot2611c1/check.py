#!/usr/bin/env python3
"""Exact checks for arXiv:2611.05000v1, Conjecture 1.

Uses only integer/rational arithmetic and exact SymPy expressions. No zeta
values, zero-free abscissa, numerical quadrature, or Lean proof are asserted.
The family check certifies the algebraic obstruction conditional on the exact
finite-additivity identity, whose analytic justification is in the draft.
"""
from fractions import Fraction
import sympy as sp


def zero(expression):
    assert sp.simplify(expression) == 0, expression


# Fidelity: source (2.1), P_{r,1}, and Theorem 2.2 reproduce source (1.15).
x, y, alpha = sp.symbols("x y alpha")
p_integrand = sum(
    sp.Integer(-1) ** (1 - m) * sp.binomial(1, m) / sp.factorial(1 + m)
    * (x - alpha) ** (1 + m) * (x - y - 1) / (x * y)
    * (1 / x + 1 / (y + 1)) ** m * (1 / (x - 1) + 1 / y)
    for m in range(2)
)
p0 = sp.simplify(sp.residue(sp.residue(p_integrand, x, 0), y, 0))
p1 = sp.simplify(sp.residue(sp.residue(p_integrand, x, 1), y, -1))
zero(p0 - alpha**2 / 2)
zero(p1 + (alpha - 1)**2 / 2)
for j in range(1, 41):
    a0 = sp.Rational(j, 4)
    value = p0.subs(alpha, a0) + (p1.subs(alpha, a0) if a0 > 1 else 0)
    expected = a0**2 / 2 if a0 < 1 else a0 - sp.Rational(1, 2)
    zero(value - expected)
print("FIDELITY: P_0,1=alpha^2/2; P_1,1=-(alpha-1)^2/2")
print("FIDELITY: source (1.15) from source (2.1) and Theorem 2.2; 40 rational alpha cases")

# Independent evaluation of the k=1 contour in Conjecture 1.
u, v = sp.symbols("u v")
z, h = sp.symbols("z h", positive=True)
f = (sp.exp(h * (u - v)) * (v - u - z) / (z * u * v)
     * (1/u + 1/(v-z) - h) * (1/(u+z) + 1/v + h))
u0 = sp.residue(f, u, 0)
um = sp.residue(f, u, -z)
zero(u0 - (z-v) * sp.exp(-h*v) / (v*z**3))
zero(um - sp.exp(-h*(v+z)) / z**2 * (1/z - 1/(v-z) + h))
main_term = sp.simplify(sp.residue(u0+um, v, 0) + sp.residue(u0+um, v, z))
zero(main_term - (1-sp.exp(-2*h*z))/z**2)
print("RESIDUES: sum=(1-exp(-2*h*z))/z^2 exactly")

# h=log(H)/2. For T=N^q, U=N^(q+1), z=2a/(qL), L=log(N),
# substitute the exact identities exp(-z log(T))=X and
# exp(-z log(U))=X*W, where X=exp(-2a), W=exp(-2a/q).
a, L, X, W = sp.symbols("a L X W", positive=True)
q = sp.symbols("q", integer=True, positive=True)
P = (1-X) * q**2 * L**2 / (4*a**2)
Q = (1-X*W) * q**2 * L**2 / (4*a**2)
D = sp.factor(Q-P)
zero(D - q**2 * X * (1-W) * L**2 / (4*a**2))
print("GAP: (Q-P)/L^2=q^2*X*(1-W)/(4*a^2); positive when 0<W<1, X>0")

# Rational Taylor enclosures: for x>0, the odd Taylor polynomial for exp(-x)
# is a lower bound and the even Taylor polynomial an upper bound, by the
# integral/Lagrange remainder. All operations below are Fraction operations.
def exp_negative_bounds(t):
    assert t > 0
    term = Fraction(1)
    upper = term
    for n in range(1, 81):
        term *= -t / n
        upper += term
    lower = upper + term * (-t) / 81
    assert lower < upper
    return lower, upper

checked_batches = 0
checked_instances = 0
positive_gap_cases = 0
for qi in range(2, 13):
    for ai in (Fraction(3, 5), Fraction(3, 4), Fraction(9, 10)):
        assert Fraction(1, 2) < ai < 1
        assert 0 < 2*ai < 3
        assert 0 < 2*ai*Fraction(qi+1, qi) < 3
        assert ai/Fraction(qi) < Fraction(1, 2)
        e_short_lower, _ = exp_negative_bounds(2*ai)
        _, e_long_upper = exp_negative_bounds(2*ai*Fraction(qi+1, qi))
        d_lower = Fraction(qi**2, 4*ai**2)*(e_short_lower-e_long_upper)
        assert d_lower > 0
        positive_gap_cases += 1
        for n in range(3, 21):
            T, U = n**qi, n**(qi+1)
            assert U == n*T
            # The chosen exponent is gamma=(q+1)/q; no real-power approximation.
            for r in range(n):
                assert 0 <= r*T < U
                checked_instances += 1
            checked_batches += 1
print(f"ADMISSIBILITY: {checked_batches} exact integer batches, {checked_instances} translated instances")
print(f"GAP_POSITIVITY: {positive_gap_cases} positive rational lower enclosures (no zeta-line assertion)")

# Exact universal linear identities for each batch size; no sample means are
# assumed to be zeta moments. Triangle inequality has coefficient sum 2,
# so |D| <= 2 max(|E_long|, |E_0|,...,|E_(n-1)|).
p, d = sp.symbols("p d", real=True)
for n in range(2, 31):
    means = sp.symbols(f"A0:{n}", real=True)
    B = sum(means) / sp.Integer(n)
    errors = [value-p for value in means]
    long_error = B-(p+d)
    zero(sum(errors)/sp.Integer(n) - long_error - d)
    coefficients = [Fraction(1, n)]*n + [Fraction(-1)]
    assert sum(abs(c) for c in coefficients) == 2
    midpoint = p+d/2
    assert all(sp.simplify(e.subs(dict.fromkeys(means, midpoint))-d/2) == 0
               for e in errors)
    zero(long_error.subs(dict.fromkeys(means, midpoint))+d/2)
print("PARTITION: exact identity and sharp factor 1/2 for batch sizes 2..30")
xi = sp.symbols("xi", positive=True)
zero(2-(1+xi)-(1-xi))
print("ASYMPTOTIC EXPONENT: 2-(1+xi)=1-xi>0 for the source range 0<xi<c<=1")
print("SCOPE: exact algebra only; analytic bridges and full family remain to be kernel-checked")
print("ALL_CHECKS_OK")
