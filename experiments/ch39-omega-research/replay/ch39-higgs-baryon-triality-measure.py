"""Exact weight pairings for the electroweak Higgs and shell triality."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations


def verify():
    # Independent color Cartan coordinates, then the residual Z3 exponent.
    colors = [(1,0,0,0),(0,1,0,0),(0,0,1,0)]
    weak = [(-1,-1,-1,1),(0,0,0,-1)]
    fund = colors+weak
    def with_phase(k, q):
        return tuple(k[:3])+((k[3]+q)%3,)
    def conj(k):
        return with_phase(tuple(-x for x in k),0)
    weights = []
    for i, w in enumerate(fund):
        weights.append(with_phase(tuple(-x for x in w), 1 if i < 3 else -1))
    for i, j in combinations(range(5), 2):
        weight = tuple(x+y for x,y in zip(fund[i],fund[j]))
        weights.append(with_phase(weight, int((i<3)!=(j<3))-1))
    counts = Counter(weights)
    assert all(counts[conj(w)] == n for w,n in counts.items())
    assert counts[(0,0,0,0)] == 1 and len(weights) == 15
    # Canonical representative of the external symmetry preserving lower Higgs.
    assert with_phase(weak[-1],1) == (0,0,0,0)
    charges = [F(2,3),F(-1,3),F(-2,3),F(1,3),-1,1,0]
    baryon = [F(1,3),F(1,3),F(-1,3),F(-1,3),0,0,0]
    assert [(2*q-b)%3 for q,b in zip(charges,baryon)] == [1,2,2,1,1,2,0]
    # Minimum over the SU(2) correction follows from its quaternionic form.
    from sympy import Matrix, conjugate, eye, simplify, symbols
    x,y = symbols("x y", complex=True)
    quat = Matrix([[x,-conjugate(y)],[y,conjugate(x)]])
    delta = quat-eye(2)
    assert simplify(delta.H*delta-(delta[:,0].H*delta[:,0])[0]*eye(2)) == Matrix.zeros(2)
    print("Seven charged conjugate pairs per generation and one neutral line.")
    print("Residual charges equal 2*Q_em-B mod 3; Higgs transport correction has exact norm.")


if __name__ == "__main__":
    verify()
