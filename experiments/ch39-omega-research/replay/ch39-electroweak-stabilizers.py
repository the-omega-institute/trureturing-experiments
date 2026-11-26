"""Weak-center parity and real structures after a parallel fundamental Higgs reduction."""
from itertools import combinations,product

import sympy as sp


def parity():
    for p in (3,5,7,11,13):
        for half in product((0,1),repeat=p//2-1):
            g = (0,1,*half)
            G = lambda k:g[min(k%p,-k%p)]
            A = lambda k:(G(k+1)+G(k-1))%2
            defects = [(A(a)+A(1)+A(a+1))%2 for a in range(p)]
            assert any(defects)
            for a in range(p):
                colors=(a,1,-a-1)
                direct=sum(G(t+1)+G(t-1) for t in colors)
                assert direct%2==defects[a]


def real_structures():
    I=sp.eye(2)
    X=sp.Matrix([[0,1],[1,0]])
    Y=sp.Matrix([[0,-sp.I],[sp.I,0]])
    Z=sp.diag(1,-1)
    kron=sp.kronecker_product
    gamma=(kron(X,I),kron(Y,I),kron(Z,X),kron(Z,Y))
    C=gamma[1]*gamma[3]
    assert C*C.conjugate()==-sp.eye(4)
    assert all(C*g.conjugate()==g*C for g in gamma)
    axes=tuple(combinations(range(4),2))
    star=sp.Matrix([[sp.LeviCivita(*a,*b) for b in axes] for a in axes])
    assert star*star==sp.eye(6)
    def wedge(U):
        return sp.Matrix([[U.extract(a,b).det() for b in axes] for a in axes])
    generators=[]
    for j in range(3):
        for u in (sp.Rational(4,5),sp.I*sp.Rational(4,5)):
            U=sp.eye(4)
            U[j,j]=U[j+1,j+1]=sp.Rational(3,5)
            U[j,j+1],U[j+1,j]=u,-sp.conjugate(u)
            generators.append(U)
    for U in generators:
        V=wedge(U)
        assert U.det()==1 and star*V.conjugate()==V*star
    assert kron(C,star)*kron(C,star).conjugate()==-sp.eye(24)


if __name__ == "__main__":
    parity()
    real_structures()
    print("Odd weak-center charge and quaternionic Weyl pairing on the fundamental-Higgs reduction verified.")
