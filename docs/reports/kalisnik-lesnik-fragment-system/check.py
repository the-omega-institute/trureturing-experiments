#!/usr/bin/env python3
"""Exact checks for Kališnik--Lešnik arXiv:2609.26560v1, Question 8.1.

No floating-point values or numerical eigenvalues are used. For a rational
configuration choose M clearing every distance denominator and set t=exp(-1/M).
The matrix has entries in Z[t]. Determinants are certified positive on 0<t<1
by factoring off endpoint roots and checking the remaining Bernstein coefficients.
This checks the stated finite cases; it is not a proof of the quantified question.
"""
from itertools import product, combinations
from math import lcm
import sympy as S

Q = S.Rational
T = S.Symbol('t')


def directions(n):
    return [s for s in product((-1, 0, 1), repeat=n) if any(s)]


def sign(x):
    return int(bool(x > 0)) - int(bool(x < 0))


def dotted(p, s, q):
    return tuple(sj if pj == qj else sign(pj-qj)
                 for pj, sj, qj in zip(p, s, q))


def center(p, r, s):
    return tuple(pj+r*sj for pj, sj in zip(p, s))


def face_distance(p, s, r, y):
    return sum(abs(yj-(pj+r*sj)) if sj else max(abs(yj-pj)-r, 0)
               for pj, sj, yj in zip(p, s, y))


def half_gap(f):
    differences = [abs(pj-qj) for p in f for q in f for pj, qj in zip(p, q) if pj != qj]
    return min(differences)/2 if differences else S.oo


def system(f, r):
    n = len(f[0]); dr = directions(n); idx = list(product(f, dr))
    distances = [[face_distance(p, s, r, center(q, r, u))
                  if dotted(p, s, q) == u else None for p, s in idx] for q, u in idx]
    m = lcm(*(S.denom(d) for row in distances for d in row if d is not None))
    matrix = S.Matrix([[0 if d is None else T**int(m*d) for d in row] for row in distances])
    return idx, m, matrix


def determinant(matrix):
    # Eliminate rows consisting of their diagonal 1. Such a row forces its
    # homogeneous coefficient to zero; the complementary principal minor has
    # exactly the original determinant, with no permutation sign.
    active = list(range(matrix.rows))
    while True:
        trivial = [i for i in active if matrix[i, i] == 1 and
                   all(matrix[i, j] == 0 for j in active if i != j)]
        if not trivial:
            break
        active = [i for i in active if i not in trivial]
    return S.Poly(matrix.extract(active, active).det(method='domain-ge'), T)


def positive_on_unit_interval(poly):
    assert not poly.is_zero
    # Extract t^a(1-t)^b, each positive on (0,1).
    for factor in (S.Poly(T, T), S.Poly(1-T, T)):
        while poly.degree() > 0 and poly.rem(factor).is_zero:
            poly = poly.exquo(factor)
    d = poly.degree()
    bernstein = [sum(poly.nth(j)*S.binomial(k, j)/S.binomial(d, j)
                     for j in range(k+1)) for k in range(d+1)]
    assert all(b >= 0 for b in bernstein), (poly.as_expr(), bernstein)
    assert any(b > 0 for b in bernstein)


def fidelity_example_7_2():
    # Source's first Section-7 example: F={(0,0),(a,0)}, a>2r.
    f = ((Q(0), Q(0)), (Q(2), Q(0)))
    for r in (Q(0), Q(1,4), Q(1,2), Q(3,4)):
        idx, m, matrix = system(f, r)
        e = T**int(m*(2-2*r))
        x = S.Matrix([e/(1+e) if (p == f[0] and s == (1,0)) or
                      (p == f[1] and s == (-1,0)) else 0 for p, s in idx])
        b = S.Matrix([sum(T**int(m*face_distance(p, (0,0), r, center(q,r,u)))
                          for p in f if tuple(sign(pj-qj) for pj,qj in zip(p,q)) == u)
                      for q,u in idx])
        assert all(S.cancel(z) == 0 for z in matrix*x-b)
    print('FIDELITY: Section 7 first example, alpha=0 and beta=e/(1+e), four rational radii: OK')


def interval_weight_identity():
    x, a, b, y = S.symbols('x a b y', real=True)
    # Three order regimes. The weight is (delta_a + Lebesgue_[a,b] + delta_b)/2.
    integrals = [S.integrate(S.exp(y-x), (x,a,b)),
                 S.integrate(S.exp(x-y), (x,a,y)) + S.integrate(S.exp(y-x), (x,y,b)),
                 S.integrate(S.exp(x-y), (x,a,b))]
    endpoint_terms = [S.exp(y-a)+S.exp(y-b), S.exp(a-y)+S.exp(y-b),
                      S.exp(a-y)+S.exp(b-y)]
    desired = [S.exp(y-a), S.Integer(1), S.exp(b-y)]
    for val, ends, target in zip(integrals, endpoint_terms, desired):
        assert S.simplify((val+ends)/2-target) == 0
    print('ANALYTIC IDENTITY: exact interval weight potential in all three order regimes: OK')


def projection_identity(f, r):
    n = len(f[0]); dr = directions(n)
    # This is the distance-addition identity underlying the homogeneous-potential
    # reduction, checked at every rational grid point in each cube below.
    sample = (Q(-1), Q(-1,2), Q(0), Q(1,2), Q(1))
    count = 0
    for q in f:
        for p,s in product(f,dr):
            u = dotted(p,s,q)
            assert any(u) and u in dr
            for v in product(sample, repeat=n):
                y = center(q,r,v)
                assert face_distance(p,s,r,y) == (
                    face_distance(p,s,r,center(q,r,u)) + face_distance(q,u,r,y))
                count += 1
    return count


def small_cases():
    cases = []
    for n, ambient in [(1, tuple(product((Q(0),Q(2),Q(4)), repeat=1))),
                       (2, tuple(product((Q(0),Q(2)), repeat=2)))]:
        for k in range(1,len(ambient)+1):
            cases.extend(tuple(c) for c in combinations(ambient,k))
    cases += [((Q(0),Q(0)),(Q(3),Q(1))),
              ((Q(0),Q(0)),(Q(2),Q(0)),(Q(-2),Q(0)),(Q(0),Q(2)),(Q(0),Q(-2))),
              ((Q(0),Q(0),Q(0)),),
              ((Q(0),Q(0),Q(0)),(Q(2),Q(0),Q(0))),
              ((Q(0),Q(0),Q(0)),(Q(2),Q(2),Q(2)))]
    determinants = projections = zero = dotted_count = 0
    for f in cases:
        g = half_gap(f)
        # Compare the coordinatewise implementation to Definition 4.5's literal
        # sign(p + chg(F)/2*s - q), using its positive-radius equivalent at infinity.
        h = Q(1) if g == S.oo else g/2
        for q, p, s in product(f, f, directions(len(f[0]))):
            literal = tuple(sign(pj+h*sj-qj) for pj,sj,qj in zip(p,s,q))
            assert dotted(p,s,q) == literal
            dotted_count += 1
        radii = (Q(0),Q(1,4),Q(1,2),Q(3,4)) if g == S.oo else (
            Q(0),g/4,g/2,3*g/4)
        for r in radii:
            assert 0 <= r < g
            idx,m,matrix = system(f,r)
            positive_on_unit_interval(determinant(matrix))
            if r == 0:
                zero += 1
            determinants += 1
            projections += projection_identity(f,r)
        print('CASE: N=%d, F=%s, four radii, determinant positive for all 0<t<1: OK'
              % (len(f[0]),f), flush=True)
    # N=0 gives Dir_0 empty and determinant of the empty matrix equal to 1.
    assert directions(0) == [] and S.det(S.zeros(0)) == 1
    print('FIDELITY: Definition 4.5 literal dotted-fragment comparison: %d indices OK' % dotted_count)
    print('FIDELITY: Lemma 5.2 (r=0) on %d configurations: OK' % zero)
    print('SETTLEMENT CHECK: %d exact determinants; %d distance identities; N=0 empty matrix: OK'
          % (determinants, projections))


if __name__ == '__main__':
    fidelity_example_7_2()
    interval_weight_identity()
    small_cases()
    print('ALL_CHECKS_OK')
