#!/usr/bin/env python3
"""Exact certificates for the Johnson--Paparella Question 9.1 candidate.

Requires Python 3 and SymPy. No approximate arithmetic is used.
Finite checks are not a Lean proof of Question 9.1 or its negation.
"""
from fractions import Fraction as F
from itertools import product
import sympy as sp

A = F(1, 4)
U = (F(1), F(1))
AXES = ((A, F(0)), (F(0), A), (-A, F(0)), (F(0), -A))


def vertex(t):
    return (-A * (1 - t*t) / (1 + t*t), -2*A*t / (1 + t*t))


def dot(p, q):
    return sum((x*y for x, y in zip(p, q)), F(0))


def cross(o, p, q):
    return (p[0]-o[0])*(q[1]-o[1]) - (p[1]-o[1])*(q[0]-o[0])


def hull(points):
    """Andrew monotone chain, with exact rational orientation tests."""
    points = sorted(set(points))
    chains = []
    for order in (points, list(reversed(points))):
        chain = []
        for p in order:
            while len(chain) >= 2 and cross(chain[-2], chain[-1], p) <= 0:
                chain.pop()
            chain.append(p)
        chains.append(chain[:-1])
    return chains[0] + chains[1]


def diamond_certificate(p):
    """Return exact convex coefficients on +x,+y,-x,-y axis vertices."""
    x, y = p
    slack = 1 - (abs(x) + abs(y))/A
    assert slack >= 0
    weights = (max(x, F(0))/A + slack/2,
               max(y, F(0))/A,
               max(-x, F(0))/A + slack/2,
               max(-y, F(0))/A)
    assert sum(weights) == 1 and all(w >= 0 for w in weights)
    assert tuple(sum((w*q[k] for w, q in zip(weights, AXES)), F(0))
                 for k in range(2)) == p
    return weights


def matrix(p):
    x, y = p
    e, f, g = sp.Matrix([1, 1, 1]), sp.Matrix([1, -1, 0]), sp.Matrix([1, 1, -2])
    return e*e.T/3 + sp.sympify(x)*f*f.T/2 + sp.sympify(y)*g*g.T/6


def fidelity():
    # Source Section 6: Walsh H_4 is ideal; C(H_4) is the row cone.
    h = sp.Matrix([[1, 1, 1, 1], [1, -1, 1, -1],
                   [1, 1, -1, -1], [1, -1, -1, 1]])
    z = sp.Matrix(sp.symbols('z0:4', real=True))
    assert h*h == 4*sp.eye(4)
    coeff = h*z/4
    mz = h*sp.diag(*z)*h/4
    for i in range(4):
        for j in range(4):
            assert sp.expand(mz[i, j] - coeff[i ^ j]) == 0
    assert h.T*coeff == z
    assert mz*sp.ones(4, 1) == sp.ones(4, 1)*z[0]
    # Thus entrywise nonnegativity is precisely coeff >= 0; row-cone
    # membership is precisely the same condition, and z0 is sum(coeff).
    assert sp.expand(sum(coeff) - z[0]) == 0
    print('FIDELITY: Walsh H4 row-cone and row-polytope descriptions reproduced exactly')


def symbolic():
    a, s, t, w = sp.symbols('a s t w', real=True)
    v = lambda r: sp.Matrix([-a*(1-r*r)/(1+r*r), -2*a*r/(1+r*r)])
    n = sp.Matrix([-(1-s*t), -(s+t)])
    c = a*(1+s*t)
    assert sp.factor(v(w).dot(v(w)) - a*a) == 0
    assert sp.factor(n.dot(v(w))-c + 2*a*(w-s)*(w-t)/(1+w*w)) == 0
    assert sp.factor(v(w)[1]/(v(w)[0]-a)-w) == 0
    alpha, beta, r0, r1, q0, q1 = sp.symbols('alpha beta r0 r1 q0 q1')
    for r, q in ((r0, q0), (r1, q1)):
        assert sp.expand((alpha+(1-alpha)*r)*(beta+(1-beta)*q) -
                         (alpha*beta + alpha*(1-beta)*q +
                          (1-alpha)*beta*r + (1-alpha)*(1-beta)*r*q)) == 0
    assert sp.expand(alpha*beta + alpha*(1-beta) + (1-alpha)*beta +
                     (1-alpha)*(1-beta) - 1) == 0
    r0, r1, q0, q1 = sp.symbols('r0 r1 q0 q1', real=True)
    assert sp.expand((r0*r0+r1*r1)*(q0*q0+q1*q1) -
                     (r0*q0+r1*q1)**2 - (r0*q1-r1*q0)**2) == 0
    print('SYMBOLIC: radius, support slack, parameter recovery, convex-product, and Lagrange identities exact')
    x, y, z = sp.symbols('x y z', real=True)
    mat = matrix((x, y))
    basis = sp.Matrix([[1, 1, 1], [1, -1, 1], [1, 0, -2]])
    assert basis.det() == 6
    assert sp.simplify(mat*basis - basis*sp.diag(1, x, y)) == sp.zeros(3)
    assert mat == mat.T
    assert mat*sp.ones(3, 1) == sp.ones(3, 1)
    charpoly = mat.charpoly()
    lam = charpoly.gen
    assert sp.factor(charpoly.as_expr() - (lam-1)*(lam-x)*(lam-y)) == 0
    assert matrix(U) == sp.eye(3)
    minimum = min(entry for p in product((-A, A), repeat=2) for entry in matrix(p))
    assert minimum == sp.Rational(1, 6)
    print('SPECTRAL: common rational eigenbasis det=6; symmetry, row sums, characteristic polynomial exact')
    print('SPECTRAL: all square-corner entries >= 1/6; A(1,1)=I')


def family(m):
    samples = [vertex(F(j, m)) for j in range(m+1)]
    base = samples + list(AXES[:2])
    points = base + [U]
    assert set(AXES) <= set(base)
    assert all(dot(p, p) == A*A for p in base)
    assert len(set(samples)) == m+1
    polygon = hull(points)
    assert len(polygon) == m+4
    edge_set = {frozenset((polygon[i], polygon[(i+1) % len(polygon)]))
                for i in range(len(polygon))}
    normals = set()
    for j in range(m):
        s, t = F(j, m), F(j+1, m)
        n, c = (-(1-s*t), -(s+t)), A*(1+s*t)
        slacks = [dot(n, p)-c for p in points]
        assert all(d <= 0 for d in slacks)
        assert {k for k, d in enumerate(slacks) if d == 0} == {j, j+1}
        assert frozenset((samples[j], samples[j+1])) in edge_set
        normal3 = (-c, *n)
        pivot = next(v for v in normal3 if v)
        normals.add(tuple(v/pivot for v in normal3))
        # Independent cone-face rank check using a nonzero 2x2 minor.
        p, q = samples[j], samples[j+1]
        assert p != q and any(q[k]-p[k] != 0 for k in range(2))
    assert len(normals) == m > 9
    pairs = 0
    for p, q in product(base, repeat=2):
        diamond_certificate((p[0]*q[0], p[1]*q[1]))
        pairs += 1
    # Convex bilinearity extends these certificates to R_m; adjoining
    # U preserves closure, and nonnegative scaling gives cone closure.
    # All generating spectra have the stated common stochastic realization.
    for p in points:
        mat = matrix(p)
        assert min(mat) >= 0
        if p != U:
            assert min(mat) >= sp.Rational(1, 6)
    return len(polygon), pairs


def main():
    fidelity()
    symbolic()
    count = total_pairs = 0
    for m in range(10, 65):
        edges, pairs = family(m)
        total_pairs += pairs
        count += 1
        if m in (10, 11, 16, 32, 64):
            print(f'FAMILY m={m}: hull_edges={edges}, nominated_facets={m}, exact_product_certificates={pairs}')
    print(f'FAMILY: {count} values m=10..64; {total_pairs} exact generator-product certificates')
    print('OBSTRUCTION INPUT: each tested polygon/cone has >9 nominated facets; a 3x3 realization has 9 entry forms')
    print('SCOPE: general facet-count argument and extension to every integer m>=10 are prose proof obligations, not Lean-checked here')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
