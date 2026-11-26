"""Closed five-link Wilson bounds and finite-volume spectral windows in the RH face model."""
from itertools import combinations
import sympy as sp

Q = sp.Rational
MASS = Q(139, 40)


def exterior_square(matrix):
    pairs = list(combinations(range(matrix.rows), 2))
    return sp.Matrix([[matrix.extract(p, q).det() for q in pairs] for p in pairs])


def involutions():
    # Symbolic two-dimensional reflection pencil, with offdiag^2=1-t^2.
    t, off, lam = sp.symbols("t off lam", real=True)
    a = sp.diag(1, -1)
    b = sp.Matrix([[t, off], [off, -t]])
    relation = {off**2: 1-t*t}
    assert (b*b).applyfunc(lambda x: sp.expand(x).subs(relation)) == sp.eye(2)
    characteristic = (lam*sp.eye(2)-a-b).det().expand().subs(relation).expand()
    assert characteristic == lam**2-2-2*t
    # Exact Gaussian-free specialization at the critical Wilson value.
    c = (4-MASS)/2
    cosine = c*c/2-1
    reflect = b.subs({t: cosine, off: sp.sqrt(1-cosine*cosine)})
    a5, b5 = [sp.diag(matrix, -1, 1, 1) for matrix in (a, reflect)]
    assert a5.det() == b5.det() == 1 and a5*a5*b5*b5 == sp.eye(5)
    for a_rep, b_rep, count in [(a5, b5, 1), (exterior_square(a5), exterior_square(b5), 3)]:
        identity = sp.eye(a_rep.rows)
        critical = (4-MASS)*identity-2*a_rep-2*b_rep
        assert len(critical.nullspace()) == count


def spectral_endpoints():
    x = sp.symbols("x", real=True)
    c = (1+sp.sqrt(5))/4
    maximal_wilson = 5*(1+c)
    real_fifth = 8*x**4-8*x*x+1
    imag_factor = 1+4*x*(2*x*x-1)
    radius2 = (5-4*x-real_fifth-MASS)**2 + Q(4, 5)*(1-x*x)*imag_factor**2
    assert radius2.subs(x, 1) == MASS*MASS
    assert sp.simplify(radius2.subs(x, -c) - (maximal_wilson-MASS)**2) == 0
    twist = sp.zeros(5)
    for i in range(4):
        twist[i, i+1] = 1
    twist[4, 0] = -1
    eigen = (2*sp.eye(5)-twist-twist.T).eigenvals()
    assert min(eigen) == 2*(1-c)
    assert 8-MASS >= MASS and maximal_wilson > 8
    a = 5-MASS
    closed_upper2 = a*a+10*a*c+25
    assert (maximal_wilson-MASS)**2 < closed_upper2 < (10-MASS)**2


def packing_constants():
    assert Q(19, 10)*(8*Q(7, 4)/56+Q(16, 64)) < 1
    print("Two-dimensional reflections supply four matter vectors, or sixteen spin-color vectors at zero.")


if __name__ == "__main__":
    involutions()
    spectral_endpoints()
    packing_constants()
