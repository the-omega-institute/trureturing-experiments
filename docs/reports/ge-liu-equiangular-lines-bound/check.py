#!/usr/bin/env python3
"""Exact checks for Ge--Liu 2606.29392v1, Conjecture 5.1 and Question 2.

Requires Python 3 and SymPy. This is an experiment, not a Lean proof.
The general graph spectrum is the classical generalized-quadrangle spectrum;
we check its polynomial consequences and independently build its q=2 instance.
"""
from fractions import Fraction
from itertools import combinations, product
from math import comb
import sympy as S


def source_fidelity():
    for alpha, expected in [(1 / (1 + 2*S.sqrt(3)), 169), (S.sqrt(5)-2, 135)]:
        bound = S.simplify((1-alpha**2)*(1-2*alpha**2)/(2*alpha**4))
        assert S.floor(bound) == expected
        print(f"SOURCE exceptional branch: {bound}; floor={expected}")
    for n in range(4, 1001):
        assert 2*n - ((n-1)*(2*n-3)//(n-2)) == 1
    print("SOURCE cycle excess identity: every n=4,...,1000 gives 1")


def small_geometry():
    # Over F_2 each nonzero vector is its own projective representative.
    points = list(product(range(2), repeat=4))[1:]
    def beta(x, y):
        return (x[0]*y[1]+x[1]*y[0]+x[2]*y[3]+x[3]*y[2]) % 2
    lines = set()
    for x, y in combinations(points, 2):
        if beta(x, y) == 0:
            z = tuple(a ^ b for a, b in zip(x, y))
            lines.add(tuple(sorted((x,y,z))))
    lines = sorted(lines)
    assert len(points) == len(lines) == 15
    B = S.Matrix([[int(x in line) for line in lines] for x in points])
    assert all(sum(B.row(i)) == 3 for i in range(15))
    assert all(sum(B.col(i)) == 3 for i in range(15))
    C = B*B.T-3*S.eye(15)
    assert C*C == 3*S.eye(15)-2*C+3*S.ones(15)
    A = S.zeros(15).row_join(B).col_join(B.T.row_join(S.zeros(15)))
    x = S.Symbol('x')
    expected = x**10*(x*x-4)**9*(x*x-9)
    assert S.expand(A.charpoly(x).as_expr()-expected) == 0
    H = (4*S.eye(30)+S.ones(30)-2*A)/5
    assert all(H[i,i] == 1 for i in range(30))
    assert all(abs(H[i,j]) == S.Rational(1,5)
               for i in range(30) for j in range(i))
    spectrum = H.eigenvals()
    assert spectrum == {S.Rational(28,5):1, S.Integer(2):1,
                        S.Rational(8,5):9, S.Rational(4,5):10, S.Integer(0):9}
    assert all(v >= 0 for v in spectrum)
    assert sum(mult for eig, mult in spectrum.items() if eig != 0) == 21
    print("GEOMETRY F_2: points=15, lines=15, degree=3, exact incidence identity OK")
    print("GEOMETRY adjacency charpoly = x^10 (x^2-4)^9 (x^2-9)")
    print("GEOMETRY Gram: 30 unit vectors, rank=21, off-diagonal magnitude=1/5; PSD by exact spectrum")


def family():
    k, q = S.symbols('k q', integer=True, positive=True)
    v = (q+1)*(q*q+1)
    m = q*(q+1)**2/2
    g = q*(q*q+1)/2
    assert S.expand(1+m+g-v) == 0
    assert S.expand((q+1)**2 + 2*q*m - v*(q+1)) == 0
    assert S.expand(2+2*m+2*g-2*v) == 0
    N = 2*v
    d = N-m
    Nk = k**6/4+k**4/2+k*k+2
    dk = 3*k**6/16+k**4/4+3*k*k/4+2
    assert S.expand(N.subs(q,k*k/2)-Nk) == 0
    assert S.expand(d.subs(q,k*k/2)-dk) == 0
    E = k**5*(k-3)/16+k**3*(k-1)/4+k*(k-3)/4+1
    assert S.expand(Nk-(dk-1+(dk-2)/k)-E) == 0
    # Uniform polynomial inequalities certified by nonnegative coefficients in k-8.
    z = S.Symbol('z', nonnegative=True)
    quartic = (4*k*k+4*k)*(4*k*k+4*k-1)/2
    for expr in [11*k**4-quartic, Nk-quartic, E-k**6/32,
                 (2*k+Nk-2*(k*k/2+1))]:
        poly = S.Poly(S.expand(expr.subs(k,z+8)),z)
        assert all(c >= 0 for c in poly.all_coeffs()) and poly.eval(0) > 0
    alpha = 1/(2*k+1)
    assert S.cancel((1-alpha**2)*(1-2*alpha**2)/(2*alpha**4)-quartic) == 0
    print("FAMILY symbolic: dimensions, multiplicities, quartic-branch equality, exact excess and E>=k^6/32 OK")
    # Field orders are prime powers 2^(2t-1); all arithmetic below is integral.
    for t in range(3, 13):
        kval = 2**t
        qval = 2**(2*t-1)
        assert 2*qval == kval*kval
        nval = 2*(qval+1)*(qval*qval+1)
        mval = qval*(qval+1)**2//2
        dval = nval-mval
        assert nval == kval**6//4+kval**4//2+kval*kval+2
        assert dval == 3*kval**6//16+kval**4//4+3*kval*kval//4+2
        assert (dval-2) % kval == 0
        quart = comb(4*kval*kval+4*kval,2)
        lin = (kval+1)*(dval-1)//kval
        eval_ = Fraction(kval**5*(kval-3),16)+Fraction(kval**3*(kval-1),4)+Fraction(kval*(kval-3),4)+1
        assert nval-lin == eval_ and eval_ >= Fraction(kval**6,32)
        assert nval > max(quart,lin)
        assert 2*kval+nval-2*(qval+1) > 0
        # Remaining nontrivial eigenvalues give positive Gram eigenvalues.
        assert all(2*kval-2*a > 0 for a in [0,-kval,-qval-1])
        if t == 3:
            assert (nval,dval,quart,lin,nval-lin) == (67650,50226,41328,56503,11147)
        print(f"t={t} k={kval} q={qval} N={nval} d={dval} quartic={quart} linear={lin} excess={nval-lin}")
    print("SCOPE: full q>=32 incidence matrices were not enumerated; family existence uses the classical symplectic construction")


if __name__ == '__main__':
    source_fidelity()
    small_geometry()
    family()
    print('ALL_CHECKS_OK')
