"""Positive comparison metrics on the 24 permutation states and two arc thresholds."""
from fractions import Fraction as Q
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb,ctx
from sympy import Matrix, Rational, eye, ones, zeros

arc = run_path(str(Path(__file__).with_name("ch39-rh-internal-arc-calibration.py")))
fixed = arc["fixed"]


def regular_comparison():
    states = list(permutations(range(4)))
    sign = [(-1)**sum(s[i]>s[j] for i in range(4) for j in range(i+1,4)) for s in states]
    sign_col = Matrix(sign)
    tstd = Matrix([[parity*(4*int(s[i]==j)-1) for i in range(4) for j in range(4)]
                   for s,parity in zip(states,sign)])
    assert tstd.rank()==9
    for mismatches,odd,expected,count in [(2,True,Rational(4,3),6), (4,True,Rational(2,3),6),
                                          (4,False,Rational(4,3),3), (3,False,1,8)]:
        steps = [s for s,k in zip(states,sign) if sum(j!=i for i,j in enumerate(s))==mismatches
                 and (k<0)==odd]
        assert len(steps)==count
        kernel=zeros(24)
        for i,s in enumerate(states):
            for v in steps:
                kernel[i,states.index(tuple(s[v[j]] for j in range(4)))]+=Rational(1,count)
        lap=eye(24)-kernel
        signval = 2 if odd else 0
        assert kernel.T==kernel and lap*ones(24,1)==zeros(24,1)
        assert lap*tstd==expected*tstd and lap*sign_col==signval*sign_col
        assert signval<=3*expected


def main():
    regular_comparison()
    ctx.prec,ctx.cap=1536,3
    data=fixed["inputs"]()
    first=arc["argument_for"](arb(33)/34,data)
    second=arc["argument_for"](arb(14)/17,data)
    d=second/(3*first)
    assert 3<d<arb(27)/8
    for length,lam in [(2,arb(2)/3),(3,arb(4)/3)]:
        u=d.root(length)
        s=(u-1)/(2-lam*u)
        assert s>0
        renorm=(1+s*lam)**length
        assert ((1+2*s)**length/renorm-d).contains(0)
        neutral_ratio=3/renorm
        neutral=arc["response"](first*neutral_ratio,data)
        print(length,"factors, difference coefficient",s.str(14),"neutral/triplet squared ratio",
              neutral_ratio.str(13),"neutral shell response",neutral.str(13))
    # Generate a dyadic coefficient for the two-factor realization. Compare
    # using Gaussian inputs directly so the common Schur pivot cancels exactly.
    nodes,weights,_=fixed["shell"]["quad"](21)
    def value(t):
        return sum(w/(1+t*z) for z,w in zip(nodes,weights))
    def invert(target):
        lo,hi=arb(0),arb(256)
        for _ in range(145):
            mid=(lo+hi)/2
            v=value(mid)-target
            if v>0:lo=mid
            elif v<0:hi=mid
            else:break
        return lo.union(hi)
    t1,t2=invert(arb(33)/34),invert(arb(14)/17)
    u=(t2/(3*t1)).sqrt()
    s=(u-1)/(2-arb(2)/3*u)
    den=2**104
    t_int=int((t1*den).mid().floor().unique_fmpz())
    s_int=int((s*den).mid().floor().unique_fmpz())
    t,s=arb(t_int)/den,arb(s_int)/den
    r=3*((1+2*s)/(1+arb(2)/3*s))**2
    assert abs(value(t)-arb(33)/34)<arb(2)**(-96)
    assert abs(value(t*r)-arb(14)/17)<arb(2)**(-96)
    print("Dyadic (t,s) numerators over 2^104:",t_int,s_int)
    print("Three-outcome neutral register response",arc["response"](first*3/8,data).str(14))


if __name__=="__main__":
    main()
