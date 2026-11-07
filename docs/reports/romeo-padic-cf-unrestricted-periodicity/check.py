#!/usr/bin/env python3
"""Exact finite checks for Romeo Conjecture 16; no infinite claim is certified here."""
from fractions import Fraction as F
from itertools import product
from math import factorial

# Q(phi), phi^2 = phi + 1. Every operation is exact.
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def neg(x): return (-x[0], -x[1])
def sub(x, y): return add(x, neg(y))
def mul(x, y):
    u, v = x; s, t = y
    return (u*s + v*t, u*t + v*s + v*t)
def inv(x):
    u, v = x; norm = u*u + u*v - v*v
    assert norm != 0
    return ((u+v)/norm, -v/norm)
def rat(q): return (F(q), F(0))
def div(x,y): return mul(x,inv(y))
ONE, PHI = rat(1), (F(0),F(1))

def vp(q,p):
    q=F(q); assert q != 0
    def vi(z):
        z=abs(z); k=0
        while z%p==0: z//=p; k+=1
        return k
    return vi(q.numerator)-vi(q.denominator)

def lift_phi(p,k):
    # Root selected by phi = 4 mod 11, or phi = 5 mod 19.
    root={11:4,19:5}[p]; mod=p
    assert (root*root-root-1)%p==0 and (2*root-1)%p!=0
    for _ in range(1,k):
        root+=next(t for t in range(p) if ((root+t*mod)**2-(root+t*mod)-1)%(mod*p)==0)*mod
        mod*=p
    assert (root*root-root-1)%mod==0
    return root

def padic_v(x,p):
    u,v=x
    if v==0:return vp(u,p)
    k=max(4,-vp(v,p)+4)
    while True:
        approx=u+v*lift_phi(p,k)
        if approx and vp(approx,p)<vp(v,p)+k:return vp(approx,p)
        k*=2

def real_interval(x,width):
    lo,hi=F(1),F(2)
    while True:
        u,v=x
        low,high=(u+v*lo,u+v*hi) if v>=0 else (u+v*hi,u+v*lo)
        if high-low<width:return low,high
        mid=(lo+hi)/2
        if mid*mid-mid-1<0:lo=mid
        else:hi=mid

def sign(x):
    if x[1]==0:return (x[0]>0)-(x[0]<0)
    width=F(1)
    while True:
        lo,hi=real_interval(x,width)
        if lo>0:return 1
        if hi<0:return -1
        width/=2

def choose(beta,n,p,variant):
    if n==0 and variant==0 and p==11:return F(14,9)
    delta=F(1,n+2)
    lo,hi=real_interval(beta,delta/8)
    left,right=max(F(0),hi-delta),lo
    assert left<right
    r=max(0,-padic_v(beta,p)); N=n+1
    k=max(4,r+N-vp(beta[1],p)+4)
    approx=beta[0]+beta[1]*lift_phi(p,k)
    z=F(p**r)*approx; modulus=p**(r+N)
    assert z.denominator%p
    c=(z.numerator*pow(z.denominator,-1,modulus))%modulus
    base=F(c,p**r)
    d=int(F(8*p**N,right-left))+1+variant
    while d%p==0:d+=1
    step=F(p**N,d)
    t=(left-base)//step+1
    a=base+step*t
    assert left<a<right
    return a

def source_fidelity():
    alphabet=list(map(F,[-2,-1,0,1,2]))+[F(1,2)]
    count=0
    # Source equation (7): A_n B_(n+1) - A_(n+1) B_n = (-1)^(n+1).
    for length in range(2,6):
        for digits in product(alphabet,repeat=length):
            Aprev,A=F(1),digits[0]; Bprev,B=F(0),F(1)
            for n,a in enumerate(digits[1:]):
                An,Bn=a*A+Aprev,a*B+Bprev
                assert A*Bn-An*B==(-1)**(n+1)
                Aprev,A,Bprev,B=A,An,B,Bn
                count+=1
    print('FIDELITY source continuant determinant: exact identities checked =',count)
    # Source [overline(1)] has quadratic equation x^2-x-1 = 0.
    assert mul(PHI,PHI)==add(PHI,ONE)
    assert sign(PHI)>0 and sign(sub(PHI,ONE))>0
    print('FIDELITY source golden mean fixed-point equation: OK')

def family_check(p,variant,steps=9):
    beta=PHI; Aprev,A=F(1),F(0); Bprev,B=F(0),F(1)
    digits=[]
    for n in range(steps):
        a=choose(beta,n,p,variant)
        diff=sub(beta,rat(a))
        assert sign(diff)>0 and sign(sub(rat(F(1,n+2)),diff))>0
        assert padic_v(diff,p)>=n+1
        if n:
            assert a>n and vp(a,p)<=-n
        digits.append(a)
        next_beta=inv(diff)
        if n==0:Aprev,A,Bprev,B=F(1),a,F(0),F(1)
        else:Aprev,A,Bprev,B=A,a*A+Aprev,B,a*B+Bprev
        assert B>0
        quotient=F(A,B)
        error=sub(PHI,rat(quotient))
        denom=mul(rat(B),add(mul(rat(B),next_beta),rat(Bprev)))
        assert error==div(rat((-1)**n),denom)
        assert padic_v(error,p)>=(n+1)**2
        if n>=1:
            assert B>factorial(n)
            assert sign(sub(rat(F(1,factorial(n)**2*(n+2))),error))>0
            assert sign(add(rat(F(1,factorial(n)**2*(n+2))),error))>0
        assert vp(B,p)<=-n*(n+1)//2
        assert sign(sub(next_beta,rat(n+2)))>0
        assert padic_v(next_beta,p)<=-(n+1)
        beta=next_beta
        print(f'FAMILY p={p} variant={variant} n={n}: a>n={a>n}; v_p(B)={vp(B,p)}; v_p(error)={padic_v(error,p)}; exact error identity OK')
    print(f'PREFIX p={p} variant={variant}: {steps} digits, first={digits[0]}, final numerator bits={digits[-1].numerator.bit_length()}')
    return digits

if __name__=='__main__':
    source_fidelity()
    first=family_check(11,0)
    second=family_check(11,1)
    third=family_check(19,1)
    assert first!=second and third
    print('Finite prefixes support the general argument; they do not certify convergence or aperiodicity of an infinite sequence.')
    print('ALL_CHECKS_OK')
