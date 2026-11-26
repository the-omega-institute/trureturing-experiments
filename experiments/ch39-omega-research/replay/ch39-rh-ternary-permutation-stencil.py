"""Minimum number of amplitudes in a positive square realizing the two-arc metric."""
from fractions import Fraction as Q
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import arb,ctx
from sympy import Rational, eye

directory=Path(__file__).parent
arc=run_path(str(directory/"ch39-rh-internal-arc-calibration.py"))
perm=run_path(str(directory/"ch39-rh-permutation-positive-edges.py"))


def parity(s):
    return (-1)**sum(i>j for i,j in combinations(s,2))


def color(p,q):
    return Q(parity(p)*parity(q)*(sum(i==j for i,j in zip(p,q))-1),3)


def main():
    states=list(permutations(range(4)))
    e=states[0]
    best=Q(0)
    triples=0
    for p,q in combinations(states[1:],2):
        # A sum-zero functional has coordinates (x,y,-x-y) on (e,p,q).
        g00=2*(1-color(e,q))
        g11=2*(1-color(p,q))
        g01=1+color(e,p)-color(e,q)-color(p,q)
        v0,v1=parity(e)-parity(q),parity(p)-parity(q)
        det=g00*g11-g01**2
        assert det>0
        bound=(g11*v0**2-2*g01*v0*v1+g00*v1**2)/det
        best=max(best,bound)
        triples+=1
    assert best==6 and triples==253
    for p in states[1:]:
        # Even arbitrary complex/quaternionic coefficients on two amplitudes:
        # 3 G_color - ss* has diagonal 2 and off-diagonal magnitude <= 2.
        off=3*color(e,p)-parity(p)
        assert abs(off)<=2
    # For 4-cycles the orbit average of ||(I-R_g)^2 f||^2 is 6I-8T_F+2T_D.
    _,double,_,four=perm["class_sums"]()
    kernel=6*eye(24)-8*four+2*double
    assert kernel.is_positive_semidefinite and kernel*four==four*kernel
    assert 6-8*(-1)+2==16
    assert 6-8*Rational(1,3)+2*Rational(-1,3)==Rational(8,3)
    ctx.prec,ctx.cap=1536,3
    data=arc["fixed"]["inputs"]()
    a=arc["argument_for"](arb(33)/34,data)
    b=arc["argument_for"](arb(14)/17,data)
    d=b/(3*a)
    s=(d-1)/(16-arb(8)/3*d)
    norm=1+arb(8)/3*s
    assert s>0 and 3<d<6 and ((1+16*s)/norm-d).contains(0)
    neutral_ratio=3/norm
    print("Best unital-difference ratios for two and three amplitudes: 3 and 6.")
    print("Ternary coefficient",s.str(14),"neutral/triplet squared ratio",neutral_ratio.str(13),
          "neutral response",arc["response"](a*neutral_ratio,data).str(14))
    print("Minimal uniform diagonal addition to restore diagonal dominance",((d-3)/2).str(14))


if __name__=="__main__":
    main()
