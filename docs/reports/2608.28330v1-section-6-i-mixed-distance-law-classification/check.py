#!/usr/bin/env python3
"""Exact helper checks for the proposed mixed distance-law classification.

Fidelity: Thäle--Tuchel, arXiv:2608.28330v1, Lemma labelled
lem:distance-identity(ii). These finite checks do not prove the classification.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import sympy as S


def compositions(total, length):
    if length == 1:
        if total > 0:
            yield (total,)
        return
    for first in range(1, total - length + 2):
        for rest in compositions(total - first, length - 1):
            yield (first,) + rest


def finite_gluing_and_weights():
    cases = bounds = triangles = 0
    for total in range(2, 10):
        for n in range(2, min(5, total) + 1):
            for comp in compositions(total, n):
                p = tuple(F(x, total) for x in comp)
                points = tuple(product(range(n), range(2)))
                internal = F(1, 8)
                radii = [F(1) + F(i, n + 1) for i in range(n)]

                def distance(x, y):
                    if x == y:
                        return F(0)
                    if x[0] == y[0]:
                        return internal
                    return radii[max(x[0], y[0])]

                observed = defaultdict(F)
                for x, y in product(points, repeat=2):
                    observed[distance(x, y)] += p[x[0]] * p[y[0]] / 4
                predicted = defaultdict(F)
                for weight in p:
                    predicted[F(0)] += weight**2 / 2
                    predicted[internal] += weight**2 / 2
                cs = [2 * p[i] * sum(p[:i], F(0)) for i in range(1, n)]
                for i, c in enumerate(cs, 1):
                    predicted[radii[i]] += c
                assert dict(observed) == dict(predicted)
                assert sum(observed.values(), F(0)) == 1
                for x, y, z in product(points, repeat=3):
                    assert distance(x, z) <= distance(x, y) + distance(y, z)
                    triangles += 1
                alpha = sum((x*x for x in p), F(0))
                delta = sum((x*x for x in p[1:]), F(0))
                C = sum(cs, F(0))
                assert alpha + C == 1
                assert p[0]**2 + delta == alpha
                for i, c in enumerate(cs, 1):
                    T = sum(p[:i], F(0))
                    assert T > 0 and p[i] == c / (2*T)
                    assert (T+p[i])**2 - T**2 == c + p[i]**2
                if alpha < 4*p[0]**2:
                    assert delta <= C*max(cs)/alpha
                    bounds += 1
                cases += 1
    print(f'FIDELITY: distance-identity(ii), {cases} rational weighted examples; '
          f'{triangles} ordered triangle checks')
    print(f'WEIGHTS: square telescoping in every example; {bounds} small-weight bounds')


def explicit_local_mixed_laws():
    # X=[0,1], uniform sampling, d(t,u)=L|t-u|. Phi(t)=R+h max(0,2t-1).
    # Then Phi_*mu = (delta_R + Uniform[R,R+h])/2 and Lip(Phi)<=1.
    r, t, u = S.symbols('r t u', real=True)
    L, h, beta = S.Rational(1, 4), S.Rational(1, 8), S.Rational(1, 32)
    assert 2*h == L
    cdf = 1 - (1-r/L)**2
    density = S.diff(cdf, r)
    assert S.simplify(density - 2*(L-r)/L**2) == 0
    assert S.integrate(density, (r, 0, L)) == 1
    assert cdf.subs(r, 0) == 0 and cdf.subs(r, L) == 1
    # Independently integrate the off-diagonal strip for each distance moment.
    for k in range(1, 9):
        direct = 2*S.integrate(S.integrate((L*(t-u))**k, (u, 0, t)), (t, 0, 1))
        via_density = S.integrate(r**k*density, (r, 0, L))
        assert S.simplify(direct-via_density) == 0
        phi_moment = S.integrate((h*(2*t-1))**k, (t, S.Rational(1, 2), 1))
        assert phi_moment == h**k/(2*(k+1))
    # Entire continuous family: density is affine decreasing, so its maximum is at 0.
    assert S.diff(density, r) < 0
    maximum_scaled_density = beta*density.subs(r, 0)
    assert maximum_scaled_density == S.Rational(1, 4)
    assert maximum_scaled_density <= S.Rational(1, 2)
    print('LOCAL MIXED 1: zeta=(delta_1+Uniform[0,1])/2; '
          'beta*nu density <= 1/4 <= 1/2; beta=1/32; nu({0})=0')
    # zeta is the restriction of delta_0/4 + 3*Exp(1)/4 to [0,1].
    # exp(-r)>=1-r: difference and first derivative vanish at zero,
    # second derivative exp(-r) is strictly positive. This is an analytic certificate.
    gap = S.exp(-r) - (1-r)
    assert gap.subs(r, 0) == 0 and S.diff(gap, r).subs(r, 0) == 0
    assert S.diff(gap, r, 2) == S.exp(-r)
    positive_r = S.symbols('positive_r', nonnegative=True)
    assert S.exp(-positive_r).is_positive
    exponential_lower_bound = S.Rational(3, 4)*(1-L)
    assert exponential_lower_bound == S.Rational(9, 16)
    assert maximum_scaled_density <= exponential_lower_bound
    print('LOCAL MIXED 2: zeta=(delta_0/4+3*Exp(1)/4)|[0,1]; '
          'beta*nu density <= 1/4 <= 9/16 <= 3*exp(-r)/4 on [0,1/4]; '
          'atom at zero untouched')
    print('LOCAL IMAGE: Phi(t)=R+(1/8)*max(0,2*t-1), '
          'Phi_*mu=(delta_R+Uniform[R,R+1/8])/2; slopes 0 and 1/4')


def atomic_donor_prefixes():
    # Uniform product measure on {0,1,2,3}^N. Distance is 4^-n at first disagreement n.
    # g(x)=sum digit_j*4^-j is uniform on [0,1]. Its tail oscillation at level n
    # is 4^-(n-1); Phi=R+max(0,2*g-1)/8 is 1-Lipschitz.
    n = S.symbols('n', integer=True, positive=True)
    q = S.Rational(3, 4)*4**(-(n-1))
    assert S.summation(q, (n, 1, S.oo)) == 1
    assert S.simplify(S.Rational(1, 4)*4**(-(n-1))-4**(-n)) == 0
    for depth in range(1, 5):
        words = tuple(product(range(4), repeat=depth))
        observed = defaultdict(F)
        for x, y in product(words, repeat=2):
            first = next((j+1 for j in range(depth) if x[j] != y[j]), None)
            observed[first] += F(1, 4**(2*depth))
        assert observed[None] == F(1, 4**depth)
        for j in range(1, depth+1):
            assert observed[j] == F(3, 4**j)
    print('ATOMIC DONOR: four-ary boundary, nu=sum_(n>=1) 3/4^n delta_(4^-n); '
          'all pair coefficients checked at depths 1..4, with exact residual 4^-depth')


def coefficient_splitting():
    for mass in [F(1, 2), F(1, 3), F(1, 7), F(5, 8)]:
        for repeats in range(1, 13):
            assert sum([mass/repeats]*repeats, F(0)) == mass
    # Repetition preserves the short ordered bin inequality, including coincident bins.
    bins = [(F(k), F(k)+F(1, 4)) for k in range(1, 9) for _ in range(3)]
    for i in range(len(bins)):
        for j in range(i+1, len(bins)):
            assert bins[i][1] <= 2*bins[j][0]
    print('SPLITTING: 48 exact atom-coefficient splits; ordered repeated-bin inequalities')


if __name__ == '__main__':
    finite_gluing_and_weights()
    explicit_local_mixed_laws()
    atomic_donor_prefixes()
    coefficient_splitting()
    print('ALL_CHECKS_OK')
