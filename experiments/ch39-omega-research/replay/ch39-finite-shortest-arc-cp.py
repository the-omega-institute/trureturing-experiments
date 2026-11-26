"""Exact finite CP test for two quadratic readings of an oriented shortest arc."""
from fractions import Fraction
from itertools import permutations
from math import gcd

from flint import acb, acb_mat, arb, ctx


def data(root5, imag, coeff):
    pts = list(permutations((0,1,3,4)))
    a,b = pts[:2]
    v = dict(zip(a,(1,(root5-1)*coeff,(1-root5)*coeff,-1)))
    G,U,D = ([[0]*3 for _ in range(3)] for _ in range(3))
    for p in pts:
        s,t = (16-sum(x*y for x,y in zip(p,z)) for z in (a,b))
        h,k = s*s+imag*s*t,t*t
        for i in range(3):
            for j in range(3):
                val = v[p[i]]*v[p[j]]
                G[i][j] += val
                U[i][j] += h*val
                D[i][j] += k*val
    return G,U,D


def disc(mat):
    d,c,b,a = mat.charpoly().coeffs()
    return (b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d).real


def certify():
    ctx.prec = 192
    G,U,D = [acb_mat(v) for v in data(arb(5).sqrt(),acb(0,1),arb(1)/2)]
    Ginv = G.inv()
    H = Ginv*U*Ginv*U.conjugate()
    K = Ginv*D*Ginv*D.conjugate()
    cp = (H*K-K*H).det()
    assert cp.imag.abs_lower()>0 and cp.real.contains(0)
    disc_h,disc_k = disc(H),disc(K)
    assert disc_h > 0 and disc_k > 0
    J = abs(cp.imag)/(2*(disc_h*disc_k).sqrt())
    assert arb(".00006122") < J < arb(".00006124")
    # Independent exact nonvanishing using a finite field specialization.
    from sympy import Matrix
    p = 101
    rt = next(j for j in range(p) if j*j%p==5)
    i = 10
    half = pow(2,-1,p)
    M,S,T = [Matrix(v) for v in data(rt,i,half)]
    _,Sc,_ = [Matrix(v) for v in data(rt,-i,half)]
    P = M.inv_mod(p)
    h,k = P*S*P*Sc,P*T*P*T
    residue = int((h*k-k*h).det())%p
    assert residue!=0
    print("CP commutator determinant mod 101",residue,"|J|",J.str(14))


if __name__=="__main__":
    certify()
