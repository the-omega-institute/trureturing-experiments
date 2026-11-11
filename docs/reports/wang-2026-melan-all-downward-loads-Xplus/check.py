#!/usr/bin/env python3
"""Exact checks for Wang's Melan heavy-load question; no floating-point evidence.

Requires Python >=3.12 and SymPy 1.14. The finite-difference tests concern
finite-dimensional models only; the continuum proof is a separate obligation.
"""

import sys
from itertools import product

import sympy as S

assert sys.version_info >= (3, 12)


def zero(expr):
    assert S.cancel(S.expand(expr.rewrite(S.exp))) == 0, expr


def fidelity():
    """Reproduce (2.2), (2.3), (2.5), and a constant-load case of (2.6)."""
    x, y = S.symbols("x y", real=True)
    cases = 0
    for r, length in [(S.Integer(1), S.Integer(1)),
                      (S.Integer(1), S.Integer(2)),
                      (S.Integer(2), S.Integer(1))]:
        a = r**2
        den = r * S.sinh(r * length)
        ka_left = S.sinh(r * x) * S.sinh(r * (length-y)) / den
        ka_right = S.sinh(r * y) * S.sinh(r * (length-x)) / den
        k0_left = x * (length-y) / length
        k0_right = y * (length-x) / length
        g_left = (k0_left-ka_left)/a
        g_right = (k0_right-ka_right)/a
        for g in [g_left, g_right]:
            zero(S.diff(g, x, 4) - a*S.diff(g, x, 2))
        zero(g_left.subs(x, 0))
        zero(S.diff(g_left, x, 2).subs(x, 0))
        zero(g_right.subs(x, length))
        zero(S.diff(g_right, x, 2).subs(x, length))
        for j in range(3):
            zero((S.diff(g_right, x, j)-S.diff(g_left, x, j)).subs(x, y))
        zero((S.diff(g_right, x, 3)-S.diff(g_left, x, 3)).subs(x, y)-1)
        # Integrate the actual two pieces of (2.5), then verify the source's
        # double integral J from the denominator of (2.6)/(2.10).
        g = S.integrate(g_right, (y, 0, x)) + S.integrate(g_left, (y, x, length))
        j = length**3/(12*a)-length/a**2+2*S.tanh(r*length/2)/r**5
        zero(S.integrate(g, (x, 0, length))-j)
        # For phi=1 and c=1, the source's rank-one correction is w=g/(1+J).
        w = g/(1+j)
        zero(S.diff(w, x, 4)-a*S.diff(w, x, 2)+j/(1+j)-1)
        zero(S.integrate(w, (x, 0, length))-j/(1+j))
        cases += 1
    print(f"FIDELITY: source Green jump, four hinged traces, double integral, "
          f"and linear constant-load solution: {cases} exact cases")


def laplacian(n):
    h = S.Rational(1, n+1)
    bmat = S.zeros(n)
    for i in range(n):
        bmat[i, i] = 2/h**2
        if i:
            bmat[i, i-1] = -1/h**2
            bmat[i-1, i] = -1/h**2
    return bmat, h


def resolvent_identities():
    cases = 0
    for n, a, b in product(range(1, 7), [1, 2], [1, 3]):
        bm, _ = laplacian(n)
        one, eye = S.ones(n, 1), S.eye(n)
        v = bm.inv()*one
        rs = [(bm+(a+b*s)*eye).inv() for s in range(4)]
        assert all(x > 0 for rr in rs for x in rr)
        for s in range(4):
            for t in range(s+1, 4):
                qs, qt = rs[s]*v, rs[t]*v
                assert qs-qt == b*(t-s)*rs[s]*rs[t]*v
                gain = (t-s)*rs[s]*rs[t]*(one+a*v)
                assert t*qt-s*qs == gain
                assert all(x > 0 for x in gain)
                assert all(x > 0 for x in qs-qt)
                cases += 1
    print(f"RESOLVENT: both identities and strict coordinate signs: {cases} exact cases")


def heavy_loads():
    s = S.symbols("s", nonnegative=True)
    # Taylor bounds for exp(1), using n! >= 7! * 8^(n-7) for n>=7.
    # This is an exact rational enclosure, not a decimal evaluation.
    e_lo = sum(S.Rational(1, S.factorial(k)) for k in range(7))
    e_hi = e_lo + S.Rational(1, S.factorial(7))*S.Rational(8, 7)
    tanh_lo, tanh_hi = (e_lo-1)/(e_lo+1), (e_hi-1)/(e_hi+1)
    j_lo, j_hi = S.Rational(1, 12)-1+2*tanh_lo, S.Rational(1, 12)-1+2*tanh_hi
    k_lo = S.Rational(5, 24)-1+S.Rational(7, 4)*tanh_lo
    assert 0 < j_lo < j_hi < S.Rational(1, 12)
    assert k_lo > 0
    assert 2*(1+j_hi)**2/(1+k_lo) < 3
    print("SOURCE_LOAD_BOUND: a=b=c=L=1, exact Taylor enclosure proves (H) threshold <3")
    for n in [2, 3, 5]:
        bm, h = laplacian(n)
        one, eye = S.ones(n, 1), S.eye(n)
        v = bm.inv()*one
        r = (bm+(1+s)*eye).inv()
        q = r*v
        g = S.cancel(h*(one.T*q)[0])
        for amplitude, profile in [(100, "uniform"), (10000, "left"),
                                   (1000000, "central")]:
            p = S.Matrix([amplitude if profile == "uniform" or
                          (profile == "left" and i < (n+1)//2) or
                          (profile == "central" and i == n//2) else 0
                          for i in range(n)])
            f = S.cancel(h*(q.T*p)[0])
            balance = S.cancel(s+s*g-f)
            numerator, denominator = S.fraction(balance)
            np, dp = S.Poly(numerator, s), S.Poly(denominator, s)
            # Positive denominator on s>=0. Sturm's algorithm counts all
            # positive roots, without requiring a coefficient sign pattern.
            assert all(c > 0 for c in dp.all_coeffs())
            assert np.eval(0) < 0
            assert np.count_roots(0, S.oo) == 1
            # Exact rational bisection of that unique zero.
            lo, hi = S.Integer(0), S.cancel(f.subs(s, 0))
            assert balance.subs(s, lo) < 0 < balance.subs(s, hi)
            start_width = hi
            for _ in range(32):
                mid = (lo+hi)/2
                if balance.subs(s, mid) < 0:
                    lo = mid
                else:
                    hi = mid
            assert hi-lo == start_width/S.Integer(2)**32
            mid = (lo+hi)/2
            rm = r.subs(s, mid)
            w = bm.inv()*rm*(p-mid*one)
            residual = S.cancel(h*(one.T*w)[0]-mid)
            assert residual == -balance.subs(s, mid)
            assert (bm*(bm+(1+mid)*eye)*w+mid*one-p) == S.zeros(n, 1)
            # Wang's (H), for a=b=c=L=1, implies ||p||_1<2(1+J)^2/(1+K).
            # Exact conservative comparison: 0<J<1/12 and K>0 imply the
            # bound is <169/72<3. Every tested discrete load has mass >3.
            assert h*sum(p) > 3
            print(f"HEAVY_DISCRETE: n={n}, profile={profile}, amplitude={amplitude}, "
                  f"mass={h*sum(p)}, positive_roots=1, "
                  f"root_bracket=[{lo},{hi}], frozen_residual=0")
    print("SCOPE: exact finite-difference checks do not prove continuum existence, "
          "uniqueness, Sobolev regularity, or convergence of Wang's Picard iteration")


def main():
    print(f"Python {sys.version_info.major}.{sys.version_info.minor}; SymPy {S.__version__}")
    fidelity()
    resolvent_identities()
    heavy_loads()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
