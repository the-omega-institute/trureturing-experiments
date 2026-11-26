"""True fractional-part pivots: dilation, separated windows and missing-prime witnesses."""
from fractions import Fraction as Q
from itertools import combinations

from flint import arb,ctx
from sympy import mobius


def true_gram(a,b,logs,bound):
    partial = arb(1)/(a*b)
    for k in range(1,bound):
        x,y=k//a,k//b
        partial += arb(1)/(a*b)-(arb(x)/b+arb(y)/a)*logs[k]+arb(x*y)/(k*(k+1))
    return partial.union(partial+arb(1)/bound)


def witness(k):
    d = [0]+[int(mobius(k//j)) if k%j==0 else 0 for j in range(1,k+1)]+[0]
    return [a-b for a,b in zip(d,d[1:])]


def pairing(u,v):
    return sum(Q(3*j*j+3*j+1,3)*a*b for j,(a,b) in enumerate(zip(u,v)))


def main():
    ctx.prec=128
    n=8192
    logs=[None]+[(arb(k+1)/k).log() for k in range(1,n)]
    G = {(a,b):true_gram(a,b,logs,n) for a in range(1,4) for b in range(a,4)}
    s23 = G[2,3]-G[1,2]*G[1,3]/G[1,1]
    s22 = G[2,2]-G[1,2]**2/G[1,1]
    s33 = G[3,3]-G[1,3]**2/G[1,1]
    corr = s23/(s22*s33).sqrt()
    assert corr>arb(".56") and corr<arb(".57")
    # Bound the weighted row sum for consecutive indices with scale ratio >=2^16.
    t=Q(1,128)
    delta=2*(2*t/(1-t)+16*t/(1-t)**2)
    assert delta<Q(3,10) and delta/(1-delta)<Q(3,7)
    for k in range(2,42):
        w=witness(k)
        assert sum(w[1:])==int(mobius(k))
        for j in range(1,44):
            integral=sum(x*(Q(2*n+1,2*j)-n//j) for n,x in enumerate(w))
            assert integral==-int(k==j)
    assert pairing(witness(2),witness(2))==16
    primes = [3,5,7,11,13,17,19]
    for p,q in combinations(primes,2):
        assert pairing(witness(p),witness(q))==Q(8,3)
    for p in primes:
        assert pairing(witness(p),witness(p))==Q(8,3)+2*p*p+Q(2,3)
    print("True sibling partial correlation after their parent is removed:",corr.str(5))
    print("Weighted separated-window leakage bound < 3/7. Omitting index 2 gives squared error at least 1/16.")


if __name__=="__main__":
    main()
