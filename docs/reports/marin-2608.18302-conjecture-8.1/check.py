#!/usr/bin/env python3
"""Exact checks for Marin arXiv:2608.18302v2, Conjecture 8.1.

Only integers, SymPy symbolic algebra, and rational evaluations are evidence.
Alternants are normalized by the same nonzero 3-by-3 Fourier determinant;
Laplace expansion makes every remaining coefficient an integer.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
import sympy as S


def sign(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))


def clean(p):
    return {a: c for a, c in p.items() if c}


def mul(p, q):
    out = defaultdict(int)
    for (a, b), c in p.items():
        for (d, e), f in q.items():
            out[a + d, b + e] += c * f
    return clean(out)


def alternant(beta):
    """det(x_i**beta_j)/D for x=(1,omega,omega²,z,z^-1,w,w^-1).

    D=det(omega**(i*j))_(0<=i,j<3) != 0.  A frozen minor is
    D*sign(residues) if its three column residues are distinct, else zero.
    This uses omega²+omega+1=0 exactly; it does not approximate omega.
    """
    out = defaultdict(int)
    for frozen in combinations(range(7), 3):
        residues = tuple(beta[j] % 3 for j in frozen)
        if len(set(residues)) != 3:
            continue
        prefactor = (-1) ** (3 + sum(frozen)) * sign(residues)
        free = [j for j in range(7) if j not in frozen]
        for p in permutations(range(4)):
            a, b, c, d = [beta[free[j]] for j in p]
            out[a - b, c - d] += prefactor * sign(p)
    return clean(out)


def orbit(a, b):
    return {(s*x, t*y) for x,y in {(a,b),(b,a)} for s,t in product((-1,1),repeat=2)}


def hull(points):
    """Andrew hull; only integer cross products; collinear interior points removed."""
    pts = sorted(set(points))
    def cross(o,a,b):
        return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    def half(seq):
        out=[]
        for p in seq:
            while len(out)>=2 and cross(out[-2],out[-1],p)<=0: out.pop()
            out.append(p)
        return out
    return tuple(half(pts)[:-1]+half(pts[::-1])[:-1])


def support(points, u, v):
    return max(u*a+v*b for a,b in points)


def evaluate(p,z,w):
    return sum(c*Fraction(z)**a*Fraction(w)**b for (a,b),c in p.items())


def main():
    z,w,T,X,Y,A,B = S.symbols('z w T X Y A B')
    # Fidelity: the frozen elementary factor is (1+T)(1+omega*T)(1+omega²*T)=1+T³.
    om=S.symbols('omega')
    frozen=S.rem(S.expand((1+T)*(1+om*T)*(1+om**2*T)),om**2+om+1,om)
    assert S.expand(frozen)==1+T**3
    D=S.expand(S.Matrix([[om**(i*j) for j in range(3)] for i in range(3)]).det())
    D=S.rem(D,om**2+om+1,om)
    assert S.rem(D**2,om**2+om+1,om)==-27
    for residues in product(range(3),repeat=3):
        minor=S.Matrix([[om**(i*j) for j in residues] for i in range(3)]).det()
        expected_minor=sign(residues)*D if len(set(residues))==3 else 0
        assert S.rem(S.expand(minor-expected_minor),om**2+om+1,om)==0
    elementary=S.Poly(S.expand((1+T**3)*(1+X*T+T**2)*(1+Y*T+T**2)),T)
    es=[S.expand(elementary.nth(i)) for i in range(8)]
    assert es==[1,X+Y,X*Y+2,X+Y+1,X+Y+1,X*Y+2,X+Y,1]
    print('Elementary coefficients:', es)
    lam=(4,4,1,1,0,0,0)
    conjugate=tuple(sum(x>=j for x in lam) for j in range(1,5))
    assert conjugate==(4,2,2,2)
    eAB=[1,A,B,A+1,A+1,B,A,1]
    J=S.Matrix(4,4,lambda i,j:eAB[conjugate[i]-i+j] if 0<=conjugate[i]-i+j<=7 else 0)
    assert J==S.Matrix([[A+1,B,A,1],[A,B,A+1,A+1],[1,A,B,A+1],[0,1,A,B]])
    det=S.expand(J.det())
    target=A**4+A**3+4*A**2+3*A-5*A**2*B-3*A*B+B**3+B**2-B
    assert det==target
    print('Dual Jacobi-Trudi determinant:',det)
    laurent=S.expand(det.subs({A:z+1/z+w+1/w,B:(z+1/z)*(w+1/w)+2}))
    coeff={}
    for term in S.Add.make_args(laurent):
        powers=term.as_powers_dict(); a=int(powers.get(z,0)); b=int(powers.get(w,0))
        c=S.cancel(term/z**a/w**b); assert c.is_Integer
        coeff[a,b]=int(c)
    expected=defaultdict(int)
    for (a,b),c in [((4,0),1),((3,3),1),((3,1),2),((3,0),1),((2,2),3),((2,0),4),((1,1),6),((0,0),10)]:
        for point in orbit(a,b):expected[point]+=c
    assert coeff==dict(expected) and len(coeff)==33
    print('F_1 exact Laurent expansion:',laurent)
    print('F_1 nonzero monomials:',len(coeff))
    P=hull(coeff); assert set(P)==orbit(4,0)|orbit(3,3)
    assert support(P,1,0)==4 and support(P,1,1)==6
    # A single signed orbit hull would force dominant generator (4,2), absent from P:
    # the integral separator (3,1) has value 14 there and maximum 12 on P.
    assert support(P,3,1)==12 and 3*4+2>12
    print('F_1 hull vertices:',P)
    print('Single-orbit obstruction: maxima force (4,2); (3,1) separates it, 14 > 12')
    delta=(6,5,4,3,2,1,0); base=(10,9,5,4,2,1,0)
    denominator=alternant(delta); numerator=alternant(base)
    assert numerator==mul(coeff,denominator)
    Q=hull(denominator)
    # Source Proposition 7.1(ii): dominant denominator vertex (N-1,N-3)=(6,4).
    assert set(Q)==orbit(6,4)
    zonotope=hull({(a+u,b+v) for a,b in product((-4,4),repeat=2) for u,v in [(2,0),(-2,0),(0,2),(0,-2)]})
    assert set(Q)==set(zonotope)
    print('FIDELITY: source denominator proposition, Q = conv W(6,4):',Q)
    for s in (1,4,7):
        sb=tuple(s*x for x in base)
        ns=alternant(sb)
        scaled={(s*a,s*b):c for (a,b),c in numerator.items()}
        assert ns==scaled
        assert len(ns)>0 and all(sb[i]>sb[i+1] for i in range(6)) and sb[-1]==0
        lam=tuple(sb[i]-delta[i] for i in range(7)); assert lam[-1]==0
        assert all(lam[i]>=lam[i+1] for i in range(6))
        # Exact numeric evaluations of the literal determinant ratio and scaling identity.
        for zz,ww in [(2,3),(3,5),(Fraction(3,2),Fraction(5,2))]:
            d=evaluate(denominator,zz,ww); assert d!=0
            lhs=evaluate(ns,zz,ww)/d
            rhs=evaluate(coeff,zz**s,ww**s)*evaluate(denominator,zz**s,ww**s)/d
            assert lhs==rhs and lhs!=0
        # Check support subtraction at a chamber grid. Newton additivity gives this
        # without expanding the high-degree Schur quotient.
        for u in range(1,13):
            for v in range(u+1):
                h=support(ns,u,v)-support(Q,u,v)
                assert h==max((10*s-6)*u+(4*s-4)*v,(9*s-6)*u+(7*s-4)*v)
        first=(10*s-6,4*s-4); second=(9*s-6,7*s-4)
        assert first[0]-second[0]==s and sum(second)-sum(first)==2*s
        assert support(ns,4,1)-support(Q,4,1)==4*first[0]+first[1]
        assert support(ns,2,1)-support(Q,2,1)==2*second[0]+second[1]
        print('Family s =',s,'lambda =',lam,'normalized numerator terms =',len(ns),
              'exposed dominant vertices =',first,second,'exact evaluations: 3/3; chamber supports: 90/90')
    print('ALL_CHECKS_OK')


if __name__=='__main__':
    main()
