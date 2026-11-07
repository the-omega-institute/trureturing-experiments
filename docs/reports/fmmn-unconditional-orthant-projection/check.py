#!/usr/bin/env python3
"""Exact polygon certificates for FMMN arXiv:2607.03582v1, Conjecture 3.

Only integer, rational, and symbolic polynomial arithmetic is used. The interval
checker proves nonnegativity of every polynomial by testing endpoints and every
stationary point (degree at most two). No sampled inequality is used as a
certificate for an interval. Requires Python 3 and SymPy.
"""
from fractions import Fraction
from itertools import combinations, product
import sympy as sp

T = sp.Symbol('t', real=True)
Q = sp.Rational
INTERVAL_CHECKS = 0
G = [(2, 1, 3, 4), (4, 4, 0, 2), (4 - T, 4, 0, 2 + T)]
A = sp.Matrix([[1, -1, 0, 1], [0, 0, 1, 1]])
POS = [(-4, 0), (4, 0), (6, 2), (6, 7), (3, 7)]
BASE = [(-10, -2), (-7, -7), (-1, -7), (6, -2),
        (10, 2), (7, 7), (1, 7), (-6, 2)]
SIGNED = [(10, 2), (10, 2 + T), (7, 7), (1, 7), (-6, 2),
          (-10, -2), (-10, -2 - T), (-7, -7), (-1, -7), (6, -2)]
# a*x + b*y <= c; all eight facets of BASE.
H = [(0, 1, 7), (0, -1, 7), (1, -1, 8), (-1, 1, 8),
     (5, 3, 56), (-5, -3, 56), (5, -7, 44), (-5, 7, 44)]


def zero(p):
    return sp.expand(p) == 0


def equal_point(a, b):
    return all(zero(x - y) for x, y in zip(a, b))


def nonnegative(p):
    """Certify p >= 0 on the CLOSED interval [0,4], with exact minima."""
    global INTERVAL_CHECKS
    p = sp.Poly(sp.expand(p), T, domain=sp.QQ)
    assert p.is_zero or p.degree() <= 2, p
    sites = [sp.Integer(0), sp.Integer(4)]
    if not p.is_zero and p.degree() == 2:
        root = -p.nth(1) / (2 * p.nth(2))
        if 0 < root < 4:
            sites.append(root)
    assert all(p.eval(x) >= 0 for x in sites), (p, sites)
    INTERVAL_CHECKS += 1


def corners(g, signed):
    return list(product(*[(-x, x) if signed else (0, x) for x in g]))


def shadow(p):
    return (sp.expand(p[0] - p[1] + p[3]), sp.expand(p[2] + p[3]))


def cross(a, b, c):
    return sp.expand((b[0] - a[0]) * (c[1] - a[1]) -
                     (b[1] - a[1]) * (c[0] - a[0]))


def area(poly):
    return sp.expand(Q(1, 2) * sum(a[0]*b[1] - a[1]*b[0]
                         for a, b in zip(poly, poly[1:] + poly[:1])))


def certify_hull(poly, points):
    # Both inclusions: each proposed vertex is an actual input point; every
    # input is in the left halfplane of every cyclic polygon edge. The same
    # halfplane checks on all proposed vertices certify convexity/order.
    assert all(any(equal_point(v, p) for p in points) for v in poly)
    for a, b in zip(poly, poly[1:] + poly[:1]):
        for p in points + poly:
            nonnegative(cross(a, b, p))
    assert area(poly).subs(T, 0) > 0
    nonnegative(area(poly))


def hull(points):
    p = sorted(set(points))
    out = []
    for seq in (p, p[::-1]):
        h = []
        for v in seq:
            while len(h) > 1 and cross(h[-2], h[-1], v) <= 0:
                h.pop()
            h.append(v)
        out += h[:-1]
    return out


def evaluate(poly, t):
    return [tuple(sp.sympify(x).subs(T, t) for x in v) for v in poly]


def main():
    for g in G:
        for coordinate in g:
            nonnegative(coordinate)
    assert all(x > 0 for x in G[0])
    assert 2**4 * sp.prod(G[0]) == 384  # contained full-dimensional box
    print('Full-dimensionality: K_t contains B(2,1,3,4), volume 384.')

    gram = A * A.T
    assert gram == sp.Matrix([[3, 1], [1, 2]]) and gram.det() == 5
    projection = A.T * gram.inv() * A
    assert projection == projection.T and projection**2 == projection
    assert A * projection == A and projection.rank() == 2
    assert sp.sqrt(gram.det()) == sp.sqrt(5)
    print('Projection: LL^T=[[3,1],[1,2]], det=5; P^2=P=P^T; LP=L.')

    positive = [shadow(v) for g in G for v in corners(g, False)]
    signed = [shadow(v) for g in G for v in corners(g, True)]
    base_points = [shadow(v) for g in G[:2] for v in corners(g, True)]
    certify_hull(POS, positive)
    certify_hull(BASE, base_points)
    certify_hull(SIGNED, signed)
    assert area(POS) == Q(87, 2)
    assert area(BASE) == 174
    assert area(SIGNED) == 174 + 3*T
    print('Interval certificate, 0<=t<=4: positive area=87/2; signed area=174+3*t.')

    # Recover BASE from its independent halfspace representation.
    intersections = []
    for (a,b,c), (d,e,f) in combinations(H, 2):
        det = a*e - b*d
        if det:
            p = (Q(c*e - b*f, det), Q(a*f - c*d, det))
            if all(u*p[0] + v*p[1] <= w for u,v,w in H):
                intersections.append(p)
    assert hull(intersections) == BASE
    u, v = (10, 2 + T), (6 - 2*T, -2 - T)
    for a,b,c in H:
        nonnegative(c - a*v[0] - b*v[1])
        nonnegative(c + a*v[0] + b*v[1])
        if (a,b,c) != (5,3,56):
            nonnegative(c - a*u[0] - b*u[1])
    assert zero(5*u[0] + 3*u[1] - 56 - 3*T)
    assert area([(10,2), u, (7,7)]) == 3*T/2
    third = [(-10,-2-T), (6-2*T,-2-T), (10,2+T), (-6+2*T,2+T)]
    certify_hull(third, [shadow(v) for v in corners(G[2], True)])
    print('Facet certificate: third shadow adds exactly two opposite triangles of area 3*t/2.')

    # Fidelity: source Proposition 4, line E=span(1,-1,0,1).
    # After the common normalization 1/sqrt(3), the full interval is [-10,10]
    # and the positive interval is [-4,6]. Thus the source upper bound is attained.
    for p in signed:
        nonnegative(p[0] + 10)
        nonnegative(10 - p[0])
    for p in positive:
        nonnegative(p[0] + 4)
        nonnegative(6 - p[0])
    for endpoint, points in [(-10,signed), (10,signed), (-4,positive), (6,positive)]:
        assert any(zero(p[0] - endpoint) for p in points)
    assert 20 == 2 * 10
    print('Fidelity (Proposition 4): line lengths 20/sqrt(3), 10/sqrt(3); ratio=2.')

    ratio = sp.cancel(area(SIGNED) / area(POS))
    assert zero(ratio - (4 + 2*T/29))
    assert zero(ratio - 4*(1 + T/58))
    assert ratio.subs(T, 2) == Q(120,29) > 4
    print('Refutation: ratio=4+2*t/29; at t=2 ratio=120/29>4.')

    # Independent rational hull implementation verifies the symbolic polygons,
    # including t=0 and t=4, on a dense finite rational grid.
    grid = sorted({Fraction(k,d) for d in range(1,13) for k in range(4*d+1)})
    for t0 in grid:
        t = Q(t0.numerator, t0.denominator)
        pos_hull = hull(evaluate(positive,t))
        full_hull = hull(evaluate(signed,t))
        assert pos_hull == hull(evaluate(POS,t))
        assert full_hull == hull(evaluate(SIGNED,t))
        assert area(pos_hull) == Q(87,2)
        assert area(full_hull) == 174 + 3*t
    print('Independent rational hull grid:', len(grid), 'parameters checked exactly.')

    for n in range(2,13):
        assert zero(ratio * 2**(n-2) - 2**n*(1+T/58))
    for k in range(1,13):
        assert zero(ratio**k - 2**(2*k)*(1+T/58)**k)
        # Bernoulli remainder has nonnegative coefficients. Its general
        # induction and the resulting unbounded amplification are in README.
        rem = sp.Poly(sp.expand((1+T/58)**k - 1 - k*T/58),T)
        assert all(x >= 0 for x in rem.all_coeffs())
    print('Products: cube padding n=2..12 and block powers k=1..12 verified exactly.')
    print('Polynomial nonnegativity certificates on [0,4]:', INTERVAL_CHECKS)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
