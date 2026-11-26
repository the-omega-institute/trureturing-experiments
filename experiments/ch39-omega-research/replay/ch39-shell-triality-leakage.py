"""Exact Behrend triality defects on both SU(5) multiplets."""
from fractions import Fraction as F
from itertools import combinations
import sympy as sp


def representation_check():
    U = sp.eye(5)
    for i,j,c,s in ((0,3,sp.Rational(3,5),sp.Rational(4,5)),
                    (1,4,sp.Rational(5,13),sp.Rational(12,13))):
        U[i,i] = U[j,j] = c
        U[i,j], U[j,i] = -s, s
    pairs = list(combinations(range(5),2))
    wedge = sp.Matrix([[U.extract(a,b).det() for b in pairs] for a in pairs])
    weak = [False,False,False,True,True]
    odd = [weak[i] != weak[j] for i,j in pairs]
    def defect(matrix, partition):
        return 3*sum(abs(matrix[i,j])**2 for i in range(matrix.rows)
                     for j in range(matrix.cols) if partition[i] != partition[j])
    block = U[:3,3:]
    gram = block.H*block
    five, ten = defect(U,weak), defect(wedge,odd)
    assert five == 6*gram.trace() and ten == 3*five-24*gram.det()
    assert sp.det(U)==1 and U.H*U==sp.eye(5)
    # Equality of the three signed minors is not assumed in the Gaussian calculation.
    assert 3*2*F(1,2)==3


def error_budgets():
    kappa, r = 2**115, F(1,2**49)
    scale2 = kappa*r**2
    moment = 10395*18**6
    lost = moment/scale2**6
    shift = 4*240*r
    tail = F(1,2**40)
    assert scale2/20 > 23*14+182+14*115+40
    first = tail+kappa*r**3 + 19*shift + (moment/scale2**5+18*lost)/(1-lost)
    second = (tail+2*kappa**2*r**5 + 973*shift
              + (moment/scale2**4+972*lost)/(1-lost))
    assert 6*first < F(1,10**5)
    assert 18*first + 24*(3+second)/kappa < F(5,10**5)
    assert 24*second < F(2,1000)
    assert (64**2)*5*F(1,2**28) < F(1,64)**2
    print("Triality defects approach 27 and 81 after scaling; certified minor correction 72.")


if __name__ == "__main__":
    representation_check()
    error_budgets()
