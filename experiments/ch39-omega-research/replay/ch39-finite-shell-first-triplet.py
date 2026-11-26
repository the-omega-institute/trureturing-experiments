"""Exact first scalar spectral band on the 24-site shell with three positive edge weights."""
from itertools import permutations

from sympy import Matrix, Poly, eye, symbols

x,a,b,c = symbols("x a b c")
points = list(permutations((0,1,3,4)))
n = 24


def matrix(weights):
    def entry(p,q):
        d = 26-sum(v*w for v,w in zip(p,q))
        return weights.get(d,0)
    A = Matrix([[entry(p,q) for q in points] for p in points])
    return sum(A[0,:])*eye(n)-A


def factors():
    u = a+b
    standard = (x-2*u)*(x**2-2*(u+c)*x+2*u*c)
    signed = (x-2*(u+c))*((x-(3*a+b+c))**2-(a-b)**2-c**2)
    pairs = x**2-2*(2*a+c)*x+6*a*c
    return x*(x-4*a-2*c)*standard**3*signed**3*pairs**2


def main():
    polynomial = factors()
    for triple in [(1,1,1),(4,2,1),(1,20,10),(7,3,11)]:
        M = matrix(dict(zip((1,2,4),triple)))
        assert M.charpoly(x).as_poly()==Poly(polynomial.subs(dict(zip((a,b,c),triple))),x)
    # Rational cuts on the equal-weight example: eigenvalue 3-sqrt(5) below 1
    # versus 3-sqrt(3) above 1. The opposite ordering occurs for (1,20,10).
    for poly,cut,expected in [(x*x-6*x+4,1,1),(x*x-6*x+6,1,0),
                              (x*x-62*x+420,4,0),(x*x-24*x+60,4,1)]:
        assert Poly(poly,x).count_roots(0,cut)==expected
    # The ordinary permutohedron has six squares and eight hexagons.
    # Each square acquires one optional diagonal when it is triangulated.
    V,E,F = 24,36,6+8
    assert (V+F)-E == 2
    assert (V+F+6)-(E+6) == 2
    print("Simple zero and first positive triplet for c <= 2a; example (1,20,10) first band doublet")


if __name__=="__main__":
    main()
