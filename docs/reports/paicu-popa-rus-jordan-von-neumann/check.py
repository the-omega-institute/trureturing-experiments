#!/usr/bin/env python3
"""Exact fidelity and reconstruction checks; finite checks do not prove the continuous claim."""
from fractions import Fraction as R
from functools import lru_cache
from itertools import product
import sympy as S

x, y, z, t = S.symbols('x y z t', real=True)
HALF = S.Rational(1, 2)


def zero(expr):
    assert S.expand(expr) == 0, S.expand(expr)


def at(g, a, b):
    return g.subs({x: a, y: b}, simultaneous=True)


def defect(f):
    return S.expand(f.subs(x, x + y) + f.subs(x, x - y) - 2*f - 2*f.subs(x, y))


def cocycle(g):
    return S.expand(at(g, x+y, z) + at(g, x-y, z) + 2*g
                    - at(g, x, y+z) - at(g, x, y-z) - 2*at(g, y, z))


def triangle_integral(expr):
    return S.expand(S.integrate((1+z)*expr, (z, -1, 0))
                    + S.integrate((1-z)*expr, (z, 0, 1)))


def reconstruct_polynomial(g):
    F = triangle_integral(at(g, x, z))
    G = triangle_integral(at(g, x, y+z))
    h = S.expand(g + defect(F)/2)
    b = S.integrate(at(g, 0, z), (z, 0, 1)) - S.integrate(at(g, 0, z), (z, -1, 0))
    a = S.expand(at(g, x, 1) + at(g, x, -1) - 2*at(g, x, 0))
    u = S.expand(-b*x/2 + HALF*S.integrate((x-t)*a.subs(x, t), (t, 0, x)))
    zero(h - (G-F))
    zero(h.subs(y, 0))
    zero(S.diff(G, y, 2) - (at(g, x, y+1)+at(g, x, y-1)-2*g))
    zero(S.diff(h, y).subs(y, 0) - b)
    zero(2*S.diff(h, y, 2) - (a.subs(x, x+y)+a.subs(x, x-y)-2*a.subs(x, y)))
    zero(S.diff(u, x, 2) - a/2)
    zero(defect(u)-h)
    f = S.expand(u-F/2)
    zero(defect(f)-g)
    return f


def polynomial_checks():
    # Version of record, Example 4.1 (v1 first example).
    g = x**2*y**2
    zero(cocycle(g))
    f_source = HALF*S.integrate((x-t)*2*t**2, (t, 0, x))
    zero(f_source-x**4/12)
    zero(defect(f_source)-g)
    print('FIDELITY: Example 4.1: Q(x^4/12) = x^2*y^2; symbolic identity verified')
    # Solve all coefficient constraints, rather than assuming g was generated as Qf.
    degree = 8
    monomials = [x**i*y**j for i in range(degree+1) for j in range(degree+1-i)]
    residuals = [S.Poly(cocycle(m), x, y, z) for m in monomials]
    powers = sorted(set().union(*(set(p.monoms()) for p in residuals)))
    matrix = S.Matrix([[p.coeff_monomial(power) for p in residuals] for power in powers])
    basis = matrix.nullspace()
    assert len(basis) == 8
    for v in basis:
        g = S.expand(sum(c*m for c, m in zip(v, monomials)))
        zero(cocycle(g))
        reconstruct_polynomial(g)
    print(f'POLYNOMIALS: all degree <= {degree} cocycles: {len(monomials)} coefficients, '
          f'constraint rank {matrix.rank()}, nullity {len(basis)}; all reconstruction identities verified')
    # Ensure the checker rejects the single-integral formula without its (x-t) factor.
    bad_u = HALF*S.integrate(2*t**2, (t, 0, x))
    assert S.expand(S.diff(bad_u, x, 2)-x**2) != 0
    print('NEGATIVE CONTROL: omitting the (x-t) factor fails u_second = a/2')


def integrate_quadratic_pieces(fun, lo, hi, breaks):
    """Exact Simpson quadrature on every polynomial piece, all of degree at most two."""
    if lo == hi:
        return R(0)
    if lo > hi:
        return -integrate_quadratic_pieces(fun, hi, lo, breaks)
    cuts = sorted({lo, hi} | {v for v in breaks if lo < v < hi})
    return sum(((b-a)*(fun(a)+4*fun((a+b)/2)+fun(b))/6
                for a, b in zip(cuts, cuts[1:])), R(0))


def corner_checks():
    grid = [R(i, 2) for i in range(-6, 7)]
    triples = pairs = 0
    for c in [R(-1, 3), R(0), R(1, 3)]:
        # q is continuous but has a genuine corner; g=Qq is not assumed smooth.
        def q(v):
            return abs(v-c)

        def g(v, w):
            return q(v+w)+q(v-w)-2*q(v)-2*q(w)

        @lru_cache(None)
        def F(v):
            return integrate_quadratic_pieces(lambda w: (1-abs(w))*g(v, w),
                                             R(-1), R(1), [R(0), c-v, v-c, c])

        b = (integrate_quadratic_pieces(lambda w:g(R(0), w), R(0), R(1), [c, -c])
             - integrate_quadratic_pieces(lambda w:g(R(0), w), R(-1), R(0), [c, -c]))

        def a(v):
            return g(v, R(1))+g(v, R(-1))-2*g(v, R(0))

        @lru_cache(None)
        def candidate(v):
            u = -b*v/2 + integrate_quadratic_pieces(lambda w:(v-w)*a(w), R(0), v,
                                                    [c-1, c, c+1])/2
            return u-F(v)/2

        for v, w, r in product(grid, repeat=3):
            assert g(v+w, r)+g(v-w, r)+2*g(v, w) == g(v, w+r)+g(v, w-r)+2*g(w, r)
            triples += 1
        for v, w in product(grid, repeat=2):
            value = candidate(v+w)+candidate(v-w)-2*candidate(v)-2*candidate(w)
            assert value == g(v, w), (c, v, w, value, g(v, w))
            pairs += 1
            # Independent structural check: candidate-q must be a homogeneous quadratic.
            k = candidate(R(1))-q(R(1))
            assert candidate(v)-q(v) == k*v*v
    print(f'CORNERS: 3 exact rational shifts of abs(x-c); {triples} cocycle triples, '
          f'{pairs} reconstruction pairs; negative, zero and positive arguments included')


if __name__ == '__main__':
    polynomial_checks()
    corner_checks()
    print('ALL_CHECKS_OK')
