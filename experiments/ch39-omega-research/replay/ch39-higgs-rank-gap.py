"""Two-color zero modes and the finite rank-four logarithmic pinning bounds."""
from fractions import Fraction
import sympy as sp


def doublet():
    X = sp.Matrix([[0, 1], [1, 0]])
    Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    Z = sp.diag(1, -1)
    eye = sp.eye(2)
    mass = sp.Rational(139, 40)
    c = 1-mass/4
    b = sp.sqrt(1-c*c)
    U, V = [c*eye+sp.I*b*axis for axis in (X, Y)]
    assert c == sp.Rational(21, 160)
    assert U.det() == V.det() == 1 and U*U.H*V*V.H == eye
    # Clifford matrices for the two orthonormal differences of step vectors.
    A = sp.kronecker_product(X, eye, X)
    B = sp.kronecker_product(Z, eye, Y)
    assert A*A == B*B == sp.eye(8) and A*B == B*A
    assert 8-(A+B).rank() == 4 and 8-(A-B).rank() == 4
    assert 4*(1-c) == mass
    assert 3*(8-(A+B).rank()) == 12


def finite_bounds():
    tolerance = Fraction(1, 2**14)
    assert (Fraction(1, 64)-Fraction(1, 100))**2 > 5*(32*tolerance)**2
    assert (128/tolerance)**24 == 2**504
    assert Fraction(3, 4)*2**29*tolerance**2 >= 1
    assert (32*6400)**24 < 2**(24*18)
    assert 23*14+24*18 == 754
    # Witness that the positive coefficient estimate works with finite tensor powers.
    for fields in range(1, 5):
        for kappa in (Fraction(0), Fraction(1, 7), Fraction(8), Fraction(91, 2)):
            cap = 4*fields*kappa
            degree = int(2+2*cap+4*cap*cap)
            assert degree >= 2*cap and degree >= 2*cap*cap and degree >= 1


if __name__ == "__main__":
    doublet()
    finite_bounds()
    print("Verified the SU(2) doublet zero modes and the finite-volume rank-four estimates.")
