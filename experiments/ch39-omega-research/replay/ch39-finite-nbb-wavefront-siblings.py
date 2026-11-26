"""Pole-cancelled local wavefronts and persistent sibling singular values."""
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

nb = run_path(str(Path(__file__).with_name("ch39-finite-nbb-pivot-boundaries.py")))


def main():
    ctx.prec=96
    n=16384
    logs=[None]+[(arb(k+1)/k).log() for k in range(1,n)]
    @lru_cache(None)
    def raw(a,b):
        return nb["true_gram"](a,b,logs,n)
    def unitary_gram(a,b):
        g=gcd(a,b)
        return arb(a*b).sqrt()/g*raw(*sorted((a//g,b//g)))
    taps = [(1,arb(1)),(4,-arb(5)/2),(16,arb(1))]
    def filtered(a,b):
        return sum(c*d*unitary_gram(i*a,j*b) for i,c in taps for j,d in taps)
    for x in (Q(2),Q(1,2)):
        assert 1-Q(5,2)*x+x*x==0
    norm2=filtered(1,1)
    assert norm2>arb("2.3") and norm2<arb("2.45")
    # Logarithmic smoothing changes each filtered column in norm by < 10^-6.
    tau,T = Q(1,2**96),2**48
    assert 2*tau*tau*T+2*tau*(1+48)+Q(1,T)<Q(1,2**46)
    assert Q(9,2)*(Q(1,2**23)+2*tau)<Q(1,10**6)
    # Hence the change in each local inner product is less than 5e-6.
    error = arb(0).union(arb("5e-6")).union(-arb("5e-6"))
    w=norm2+error
    ab,ac,bc=[(filtered(i,j)+error)/w for i,j in ((1,2),(1,3),(2,3))]
    local = bc-ab*ac
    assert local>arb(".33") and local<arb(".36") and w>2
    # Absolute inter-cluster row sum, including every other cluster at all levels.
    t=Q(1,2**32)
    eps=6*102*(7*t/(1-t)+64*t/(1-t)**2)
    assert eps<Q(1,50000) and 6*eps/(1-eps)<Q(1,1000)
    print("Pole moments zero. Smoothed local residual",local.str(5),
          "; all K cross-branch singular values lie between .329 and .361 uniformly in K.")


if __name__=="__main__":
    main()
