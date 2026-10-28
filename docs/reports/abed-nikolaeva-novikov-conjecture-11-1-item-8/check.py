#!/usr/bin/env python3
"""Exact checks for arXiv:2610.34388v1, Conjecture 11.1, item 8.

Uses Python integers/Fraction and SymPy exact matrices. Finite checks are evidence
for the proposed proof route, not a proof of the universal matrix assertion.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
import sympy as S


def omc_diagonal(xs):
    assert xs and min(xs) >= 0 and max(xs) > 0
    return (max(xs) - min(xs)) / (max(xs) + min(xs))


def fidelity():
    k, ell = S.symbols("k ell", positive=True)
    rho = S.diag(k, 1)
    A = S.Matrix([[0, 1], [ell, 0]])
    absA = S.diag(ell, 1)
    absAAstar = S.diag(1, ell**2)
    assert absA**2 == A.adjoint() * A
    assert absAAstar == A * A.adjoint()
    assert absAAstar**2 == (A * A.adjoint()).adjoint() * (A * A.adjoint())
    L = S.trace(rho * absAAstar)
    R = S.trace(rho * A.adjoint() * A)
    assert L == k + ell**2 and R == 1 + k * ell**2
    d, e = S.symbols("d e", real=True)
    expr = (L/R - 1).subs({k: (1+d)/(1-d), ell: (1-e)/(1+e)})
    assert S.factor(expr - 4*d*e/(1+e**2-2*d*e)) == 0
    count = 0
    for d, e, r in product([F(0),F(1,5),F(3,5),F(9,10)],
                            [F(0),F(1,4),F(1,2),F(3,4),F(1)], range(1,5)):
        k = (1+d)/(1-d)
        ell = (1-e)/(1+e)
        assert omc_diagonal([k,1]*r) == d
        assert omc_diagonal([ell,1]*r) == e
        L, R = r*(k+ell*ell), r*(1+k*ell*ell)
        assert R > 0 and L/R-1 == 4*d*e/(1+e*e-2*d*e)
        count += 1
    k, ell = F(4), F(1,3)
    L, R = k+ell*ell, 1+k*ell*ell
    assert (L,R,L/R-1) == (F(37,9),F(13,9),F(24,13))
    print("FIDELITY: Proposition 5.22 symbolic block and parameter identity OK")
    print(f"FIDELITY: {count} repeated-block cases, n=2,4,6,8, endpoints included OK")
    print(f"EXACT_INSTANCE: delta=3/5 eta=1/2 L={L} R={R} violation={L/R-1}")


def scalar_identities():
    x,y,K,Y = S.symbols("x y K Y")
    f = lambda a,b: (a*b+1)/(a+b)
    assert S.factor(f(K,y)-f(x,y)-(K-x)*(y*y-1)/((K+y)*(x+y))) == 0
    assert S.factor(f(x,Y)-f(x,y)-(Y-y)*(x*x-1)/((x+Y)*(x+y))) == 0
    print("ORBIT_FACTOR: both scalar monotonicity identities OK")


def small_spectra():
    count=0
    for n in range(2,6):
        spectra = list(combinations_with_replacement([F(1),F(2),F(3)], n))
        for rr, hh in product(spectra, repeat=2):
            r,h=rr[::-1],hh[::-1]
            k,t=r[0]/r[-1],h[-1]/h[0]
            C=(k+t)/(1+k*t)
            lo=sum(a*b for a,b in zip(r,h[::-1]))
            hi=sum(a*b for a,b in zip(r,h))
            assert hi <= C*lo
            for hp in set(permutations(h)):
                v=sum(a*b for a,b in zip(r,hp))
                assert lo <= v <= hi and v <= C*lo
                count += 1
    print(f"ORBIT_CHECKS: {count} exact permutation pairings for all spectra in {{1,2,3}}, n=2..5 OK")


def rational_matrix_checks():
    # Rational complex unitaries permit independent left/right singular frames.
    I=S.eye(2)
    Q=S.Matrix([[S.Rational(3,5),-S.Rational(4,5)],
                [S.Rational(4,5),S.Rational(3,5)]])
    phase=S.diag(S.I,1)
    frames=[I,Q,phase*Q,Q*phase]
    assert all(U.adjoint()*U == I for U in frames)
    count=zero=0
    for di,ei,i,j in product(range(4),range(5),range(4),range(4)):
        d=[S.Rational(0),S.Rational(1,5),S.Rational(3,5),S.Rational(9,10)][di]
        e=[S.Rational(0),S.Rational(1,4),S.Rational(1,2),S.Rational(3,4),S.Rational(1)][ei]
        k,ell=(1+d)/(1-d),(1-e)/(1+e)
        T=S.diag(ell,1)
        rho=frames[(i+2*j)%4]*S.diag(k,1)*frames[(i+2*j)%4].adjoint()
        A=frames[i]*T*frames[(j+1)%4].adjoint()
        B=S.Rational(1+(i+j)%3,2)*frames[j]*T*frames[(2*i+1)%4].adjoint()
        X=A*B.adjoint()
        H=(X.adjoint()*X).applyfunc(S.expand)
        # det(H) = |det(X)|^2 has an exact nonnegative rational square root.
        det=S.expand(H.det())
        root=S.sqrt(det)
        assert root.is_Rational and root >= 0
        ss=S.expand(S.trace(H)+2*root)
        N=S.expand(S.trace(rho*(A.adjoint()*A+B.adjoint()*B))/2)
        assert N.is_Rational and N > 0
        C=(k+ell**2)/(1+k*ell**2)
        if ss == 0:
            assert H == S.zeros(2) and X == S.zeros(2)
            zero += 1
        else:
            assert ss.is_Rational and ss > 0
            Z=H+root*I
            # For H >= 0, sqrt(H)=(H+sqrt(det H) I)/sqrt(tr H+2sqrt(det H)).
            # Z is PSD: its principal minors are nonnegative real rationals.
            assert (Z**2-ss*H).applyfunc(S.expand) == S.zeros(2)
            assert Z[0,0] >= 0 and Z[1,1] >= 0 and S.expand(Z.det()) >= 0
            numerator=S.expand(S.trace(rho*Z))
            assert numerator.is_Rational and numerator >= 0
            # Squaring is equivalent here: both sides are nonnegative.
            assert numerator**2 <= (C*N)**2*ss
        count += 1
    print(f"MATRIX_CHECKS: {count} independently oriented exact rational-complex 2x2 inputs OK ({zero} zero products)")


if __name__ == "__main__":
    fidelity()
    scalar_identities()
    small_spectra()
    rational_matrix_checks()
    print("ALL_CHECKS_OK")
