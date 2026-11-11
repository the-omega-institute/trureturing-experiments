#!/usr/bin/env python3
"""Exact certificates for the Arts pointwise-necessity counterexample family.

All evidence uses integers, Fraction, and exact SymPy expressions. The analytic
steps justified in README.md (concavity, logarithmic inequalities, Poisson
moments, Chebyshev, and a union bound) are not inferred from finite sampling.
"""

from fractions import Fraction as Q
from itertools import product
import sympy as s


def weights(z):
    w = [Q(1) for _ in z]
    for a in (0, 1):
        ids = [i for i, b in enumerate(z) if b == a]
        for j in range(0, len(ids) - 1, 2):
            w[ids[j]], w[ids[j + 1]] = Q(99, 100), Q(101, 100)
    return w


def nonnegative_rational(expr, var, strict=False):
    """Certificate valid for every var >= 0: numerator/denominator coefficients."""
    num, den = s.fraction(s.cancel(expr))
    pn, pd = s.Poly(num, var), s.Poly(den, var)
    if pd.eval(0) < 0:
        pn, pd = -pn, -pd
    assert all(x >= 0 for x in pn.all_coeffs()), pn
    assert all(x >= 0 for x in pd.all_coeffs()) and pd.eval(0) > 0, pd
    if strict:
        assert pn.eval(0) > 0, pn


def fidelity():
    # Source equations (E_oab), (E_Da), including doubled Poisson self-loops.
    n, rho = 4, Q(1, 1000)
    mat = ((Q(64), Q(128)), (Q(128), Q(257)))
    cases = 0
    for z in product((0, 1), repeat=n):
        w = weights(z)
        for e in product((0, 1), repeat=n):
            r = [[sum((w[i] for i in range(n) if e[i] == a and z[i] == b), Q(0)) / n
                  for b in (0, 1)] for a in (0, 1)]
            eo = [[sum((rho * w[i] * w[j] * mat[z[i]][z[j]]
                        for i in range(n) for j in range(n)
                        if e[i] == a and e[j] == b), Q(0))
                   for b in (0, 1)] for a in (0, 1)]
            for a, b in product((0, 1), repeat=2):
                rhs = rho * n * n * sum((r[a][c] * mat[c][d] * r[b][d]
                                         for c, d in product((0, 1), repeat=2)), Q(0))
                assert eo[a][b] == rhs
            for a in (0, 1):
                ed = sum((rho * w[i] * w[j] * mat[z[i]][z[j]]
                          for i in range(n) for j in range(n) if e[i] == a), Q(0))
                assert sum(eo[a], Q(0)) == ed
            cases += 1
    print(f"FIDELITY: source E[o_ab]=rho*n^2*(R*S*R^T)_ab and E[D_a], {cases} cases")


def normalization():
    cases = 0
    for n in range(11):
        for z in product((0, 1), repeat=n):
            w = weights(z)
            assert all(Q(99, 100) <= x <= Q(101, 100) for x in w)
            for a in (0, 1):
                assert sum((w[i] for i in range(n) if z[i] == a), Q(0)) == z.count(a)
            if n >= 3:
                assert Q(99, 100) in w and Q(101, 100) in w
            cases += 1
    print(f"WEIGHTS: all {cases} binary labelings, n=0..10; normalized, common support, nonconstant for n>=3")


def family():
    u = s.Symbol('u', nonnegative=True)
    L = u + 64
    mat = s.Matrix([[L, 2*L], [2*L, 4*L+1]])
    assert s.expand(mat.det() - L) == 0
    for entry in mat:
        nonnegative_rational(entry, u, strict=True)
    alpha, beta = 3*L/2, (6*L+1)/2
    p = [s.Rational(1, 3), s.Rational(2, 3)]
    q = [2*L/(6*L+1), (4*L+1)/(6*L+1)]
    assert s.cancel(sum(p)-1) == s.cancel(sum(q)-1) == 0
    chi = sum((p[i]-q[i])**2/q[i] for i in (0, 1))
    bound = 1/(12*(4*L+1))
    assert s.cancel(alpha*chi-bound) == 0
    nonnegative_rational(1-bound, u, strict=True)
    # For r=beta/alpha>=2, log r >= 2(r-1)/(r+1) implies
    # KL(alpha||beta)-KL(beta||alpha) <= 0.
    r = beta/alpha
    nonnegative_rational(r-2, u, strict=True)
    x = (r-1)/(r+1)
    nonnegative_rational(x, u, strict=True)
    nonnegative_rational(1-x, u, strict=True)
    print("FAMILY: det(S_L)=L>0; alpha=3L/2; beta=(6L+1)/2; L=64+u, u>=0")
    print("CONSTANT: alpha*chi2(p||q)=1/[12(4L+1)]<1, symbolic for all L>=64")

    # Symbolic derivative, concavity and KL decomposition, not sampled t-values.
    t = s.Symbol('t', real=True)
    a0, a1, b0, b1 = s.symbols('a0 a1 b0 b1', positive=True)
    A, B = (a0+a1)/2, (b0+b1)/2
    kl = lambda x, y: x*s.log(x/y)+y-x
    h = lambda x, y: (1-t)*x+t*y-x**(1-t)*y**t
    f = (h(a0,b0)+h(a1,b1))/2-t*kl(A,B)
    deriv = s.diff(f,t).subs(t,0)
    normalized_kl = A*sum((aa/(2*A))*s.log((aa/(2*A))/(bb/(2*B)))
                          for aa,bb in [(a0,b0),(a1,b1)])
    assert s.simplify(s.expand_log(deriv-normalized_kl, force=True)) == 0
    curvature = -sum(aa**(1-t)*bb**t*s.log(bb/aa)**2
                     for aa,bb in [(a0,b0),(a1,b1)])/2
    assert s.simplify(s.expand_log(s.diff(f,t,2)-curvature, force=True)) == 0
    f_printed = (h(a0,b0)+h(a1,b1))/2-t*kl(B,A)
    assert s.simplify(s.expand_log(
        s.diff(f_printed,t).subs(t,0)-normalized_kl-kl(A,B)+kl(B,A), force=True)) == 0
    print("DERIVATIVES: both printed H orientations checked symbolically; curvature is -sum(positive*log^2)/2")

    # Uniform conditional degree bounds on |n1/n-1/2|<=1/100, n>=400.
    lo, hi = s.Rational(99,100), s.Rational(101,100)
    upper = hi*s.Rational(151,100)*L
    nonnegative_rational(s.Rational(25,16)*L-upper, u, strict=True)
    lower = lo*s.Rational(149,50)*L-hi**2*(4*L+1)/400
    nonnegative_rational(lower-s.Rational(47,16)*L, u, strict=True)
    # n>=400 makes the negative diagonal correction no larger than at 400.
    print("MEANS: mu_1/log n <=25L/16; mu_2/log n >=47L/16, n>=400 on the size event")
    assert Q(1,4)/Q(1,100)**2 == 2500
    log54_lower = Q(1,4)/(1+Q(1,4))
    log43_upper = Q(1,4)+Q(1,4)**2/(2*(1-Q(1,4)))
    assert log54_lower == Q(1,5)
    assert log43_upper == Q(7,24)
    e1 = Q(25,16)/4-Q(35,16)*log54_lower
    e2 = -Q(47,16)/4+Q(35,16)*log43_upper
    assert e1 == -Q(3,64) and e2 == -Q(37,384)
    nonnegative_rational(3*L/64-1, u, strict=True)
    nonnegative_rational(37*L/384-1, u, strict=True)
    print("TAILS: exact log certificates give -3L/64 and -37L/384; size-event failure <=2500/n")
    print("UNION: failure <=2500/n+n^(1-3L/64)+n^(1-37L/384), all sufficiently large n")
    assert s.cancel(bound.subs(u,0)) == s.Rational(1,3084)
    assert Q(35*64,16) == 140
    assert 1+e1*64 == -2 and 1+e2*64 == -Q(31,6)
    print("L=64: S=((64,128),(128,257)); threshold=140 log n; C<=1/3084; exponents=-2,-31/6")


if __name__ == '__main__':
    fidelity()
    normalization()
    family()
    print('ALL_CHECKS_OK')
