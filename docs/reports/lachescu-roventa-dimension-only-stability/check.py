#!/usr/bin/env python3
"""Exact checks for arXiv:2207.08058v1, Problem 1.

Finite checks and symbolic identities support the accompanying general argument;
they are not a Lean proof or an exhaustive verification over real functions.
No floating-point calculation is used.
"""
from fractions import Fraction as F
import sympy as S


def fidelity():
    # Section 2: every convex f with f(0) <= 0 is m-star-convex.
    # Reproduce the barycentric coefficients used in the source's proof.
    m, l, x, y = S.symbols('m l x y', real=True)
    assert S.expand((1-l) + m*l + (l-m*l) - 1) == 0
    assert S.expand((1-l)*x + m*l*y + (l-m*l)*0 - ((1-l)*x+m*l*y)) == 0
    # For f(t)=t^2, an exact nonnegative gap certificate on 0<=m,l<=1.
    gap = (1-l)*x**2 + m*l*y**2 - ((1-l)*x+m*l*y)**2
    cert = l*(1-l)*(x-m*y)**2 + m*l*(1-m)*y**2
    assert S.expand(gap-cert) == 0
    # Also reproduce the source's displayed nonconvex polynomial example.
    def f(t): return t**4-5*t**3+9*t**2-5*t
    nonconvex_gap = f(F(5,4)) - (f(F(1))+f(F(3,2)))/2
    assert nonconvex_gap == F(5,256) > 0
    samples = 0
    for a in range(17):
        for b in range(17):
            for c in range(1,8):
                u, v, t, mm = F(a,4), F(b,4), F(c,8), F(16,17)
                assert f((1-t)*u+mm*t*v) <= (1-t)*f(u)+mm*t*f(v)
                samples += 1
    print('FIDELITY: source Section 2 barycentric proof coefficients and square-gap identity OK')
    print(f'FIDELITY: source polynomial nonconvex midpoint gap={nonconvex_gap}; m=16/17 samples={samples} OK')


def family_checks():
    m, l, x, y, w = S.symbols('m l x y w', real=True)
    # w denotes omega(|x-y|); no assumption on w is used in these identities.
    z = (1-l)*x+m*l*y
    affine_gap = (1-l)*x+m*l*y-z-m*l*(1-l)*w
    assert S.expand(affine_gap + m*l*(1-l)*w) == 0
    shift = -m/(1-m)
    shifted_gap = (1-l)*(x+shift)+m*l*(y+shift)-(z+shift)-m*l*(1-l)*w
    assert S.cancel(shifted_gap-m*l*(1-(1-l)*w)) == 0
    # For omega(t)=sqrt(t), compare nonnegative sides by squaring, exactly.
    # lambda(1-lambda)*sqrt(d)<=1/4 gives approximate membership with delta=m/4;
    # (1-lambda)*sqrt(d)<=1 gives exact membership of x-m/(1-m).
    cases = 0
    for ni in range(1,17):
        mm = 1-F(1,2*ni)
        assert 0 < mm < 1 and 0 < mm/4 <= F(1,4)
        for a in range(9):
            for b in range(9):
                d = abs(F(a,8)-F(b,8))
                for c in range(1,8):
                    t = F(c,8)
                    zz = (1-t)*F(a,8)+mm*t*F(b,8)
                    assert 0 <= zz <= 1
                    assert 16*t*t*(1-t)**2*d <= 1
                    assert (1-t)**2*d <= 1
                    cases += 1
    print(f'FAMILY: affine and exact-comparator gap identities OK; squared radical grid cases={cases} OK')

    # The special pair x=my makes the interpolation point equal to x.
    u, v, a, b = S.symbols('u v a b', real=True)
    raw = (1-l)*u+m*l*v-m*l*(1-l)*a*b-u
    factored = l*(m*v-u-m*(1-l)*a*b)
    assert S.expand(raw-factored) == 0
    print('ONE_STEP: x=my cancellation and lambda-factor identity OK (no continuity of h invoked)')

    for ni in range(1,17):
        hs = S.symbols(f'h0:{ni+1}')
        telescoped = sum(m**(ni-1-j)*(hs[j+1]-m*hs[j]) for j in range(ni))
        assert S.expand(telescoped-(hs[ni]-m**ni*hs[0])) == 0
    print('TELESCOPE: arbitrary symbolic h-values, lengths 1..16 OK')

    checks = 0
    for ni in range(1,65):
        mm = 1-F(1,2*ni)
        assert mm**ni >= F(1,2)  # Bernoulli bound used by the general proof.
        for p,r in [(1,2),(1,3),(2,3),(3,4)]:
            assert 0 < p < r
            # All terms m^(n-j+j*q) >= m^n: compare their positive r-th powers.
            for j in range(ni):
                assert mm**(r*(ni-j)+j*p) >= mm**(r*ni)
                checks += 1
            # The weak deficit, norm and k bounds, compared by positive r-th powers.
            deficit_power = F(ni**r,2**r)*(1-mm)**p
            assert deficit_power == F(ni**(r-p),2**(r+p))
            norm_bound_power = deficit_power/2**r
            k_bound_power = 4**r*norm_bound_power
            assert k_bound_power == F(ni**(r-p),2**p)
    print(f'LOWER_BOUNDS: n=1..64; q=1/2,1/3,2/3,3/4; positive-power term comparisons={checks} OK')
    print('SQUARE_ROOT: M^2 >= n/32 and k_1^2 >= n/2 are the derived weak bounds')
    for k in range(1,65):
        ni = 2*k*k+1
        assert F(ni,2) > k*k
    print('UNBOUNDED_FAMILY: integer proposed constants k=1..64 defeated by n=2*k^2+1 exactly')


if __name__ == '__main__':
    fidelity()
    family_checks()
    print('ALL_CHECKS_OK')
