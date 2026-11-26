"""An interpolating pole-cancelled wavefront for two sparse prime subtrees."""
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

nb = run_path(str(Path(__file__).with_name("ch39-finite-nbb-pivot-boundaries.py")))


def main():
    ctx.prec=128
    n=16384
    logs=[None]+[(arb(j+1)/j).log() for j in range(1,n)]
    @lru_cache(None)
    def bare(a,b):
        return nb["true_gram"](a,b,logs,n)
    def U(a,b):
        k=gcd(a,b)
        return arb(a*b).sqrt()/k*bare(*sorted((a//k,b//k)))
    taps=[(1,arb(1)),(4,-arb(5)/2),(16,arb(1))]
    @lru_cache(None)
    def local(a,b):
        error=arb("5e-6")
        return sum(c*d*U(a*i,b*j) for i,c in taps for j,d in taps)+(-error).union(error)
    norm0=local(1,1)
    assert norm0>arb("2.3") and norm0<arb("2.45")
    assert abs(local(2,3)/norm0)<arb(".5")
    assert arb("2.3")*(1-arb(".5")**2)>arb("1.7")
    for s in (0,1):
        def gram(a,b):
            return (1+s*s)*local(a,b)-s*(local(2*a,3*b)+local(3*a,2*b))
        mass=gram(1,1)
        c=gram(2,3)/mass-gram(1,2)*gram(1,3)/mass**2
        assert (s==0 and c>arb(".32")) or (s==1 and c<-arb(".34"))
    t=Q(1,2**32)
    error=18*700*(7*t/(1-t)**2+64*t*(1+t)/(1-t)**3)
    assert error<Q(1,4000) and 6*error/(1-error)<Q(2,1000)
    print("Two smooth wavefront endpoints straddle zero local residual; uniform cross-branch trace-norm budget below .002.")


if __name__=="__main__":
    main()
