"""A finite two-layer Wilson index on the 24-point permutation shell."""
from pathlib import Path
from runpy import run_path

from sympy import Matrix, Poly, Rational as Q, eye, symbols

shell = run_path(str(Path(__file__).with_name("ch39-finite-shell-first-triplet.py")))
x = symbols("x")
L = shell["matrix"]({1:1,2:1,4:1})
I = eye(24)


def kernel(lo,hi,t):
    return (L-hi*I).row_join(t*I).col_join((t*I).row_join(lo*I-L))


def main():
    H = kernel(Q(1,2),Q(1),Q(1,8))
    P = H.charpoly(x).as_poly()
    factors = P.factor_list()[1]
    def count(left,right):
        return sum(m*f.count_roots(left,right) for f,m in factors)
    # Gershgorin for graph degree four and these mass shifts puts all roots in (-10,10).
    assert count(-10,-Q(3,25))==27 and count(Q(3,25),10)==21
    assert -Q(21-27,2)==3
    # The unique first-band polynomial is irreducible over Q.
    band = Poly(x*x-6*x+4,x)
    assert band.is_irreducible
    assert band.count_roots(Q(19,25),Q(77,100))==1
    assert (1-Q(77,100))*(Q(19,25)-Q(1,2)) > Q(1,64)
    # Exact three-dimensional target eigenspace over Q(sqrt(5)).
    assert (L*L-6*L+4*I).rank()==18
    # A same-mass model has H² positive with lower bound t² and zero index.
    H0 = kernel(Q(1),Q(1),Q(1,8))
    block = (L-I)**2+I/64
    zeros = Matrix.zeros(24)
    assert H0**2 == block.row_join(zeros).col_join(zeros.row_join(block))
    print("Hermitian Wilson spectrum below -3/25: 27, above +3/25: 21; index 3, no opposite zero")


if __name__=="__main__":
    main()
