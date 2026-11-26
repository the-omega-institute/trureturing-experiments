"""Finite quaternion heap measurements of the four-cycle second difference."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

arc=run_path(str(Path(__file__).with_name("ch39-rh-internal-arc-calibration.py")))


def norm2(q):
    return sum(v*v for v in q)


def conj(q):
    return (q[0],-q[1],-q[2],-q[3])


def mul(x,y):
    a,b,c,d=x
    w,p,q,r=y
    return (a*w-b*p-c*q-d*r,a*p+b*w+c*r-d*q,a*q-b*r+c*w+d*p,a*r+b*q-c*p+d*w)


def gap2(a,b):
    return norm2(tuple(x-y for x,y in zip(a,b)))


def samples():
    states=list(permutations(range(4)))
    parity={p:(-1)**sum(i>j for i,j in combinations(p,2)) for p in states}
    four=[p for p in states if parity[p]<0 and all(i!=j for i,j in enumerate(p))]
    h=((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))
    columns=[parity]+[{p:parity[p]*sum(h[t][i]*h[p[t]][j] for t in range(4))//4 for p in states}
                       for i in range(3) for j in range(3)]
    paths=[]
    for p,g in product(states,four):
        q=tuple(p[i] for i in g)
        paths.append((p,q,tuple(q[i] for i in g)))
    result=[]
    for col in columns:
        hist=Counter((col[p],col[q],col[r]) for p,q,r in paths)
        result.append((sum(v*v for v in col.values()),hist))
    return result


def main():
    ctx.prec,ctx.cap=1536,3
    data=arc["fixed"]["inputs"]()
    target=arc["argument_for"](arb(14)/17,data)/(3*arc["argument_for"](arb(33)/34,data))
    probes=samples()
    shell=[s for s in product(range(5),repeat=3) if sum(i*i for i in s)==17]
    rotors={}
    for a,b in combinations(shell,2):
        raw=mul(conj((0,*a)),(0,*b))
        rotors[Q(raw[0],17)]=tuple(Q(v,17) for v in raw)
    passed=[]
    for x,rotor in sorted(rotors.items()):
        unit=(Q(1),0,0,0)
        assert norm2(rotor)==1
        chords=gap2(unit,rotor)
        lookup={1:rotor,0:unit,-1:conj(rotor)}
        values=[]
        for scale,hist in probes:
            total=0
            for (i,j,k),count in hist.items():
                v,w,z=lookup[i],lookup[j],lookup[k]
                pred=mul(mul(w,conj(v)),w)
                diff=gap2(pred,z)
                assert diff==gap2(mul(conj(v),w),mul(conj(w),z))
                total+=count*diff
            values.append(total/(6*scale*chords))
        sign,color=values[:2]
        assert sign==8*x*x*(1+x) and values[1:]==[Q(4,3)*(1+x)]*9
        assert sign/color==6*x*x
        def to_arb(q):return arb(q.numerator)/q.denominator
        if to_arb(sign/color)>target:
            passed.append(x)
            s=(target-1)/(to_arb(sign)-target*to_arb(color))
            assert s>0
            print("Arc", (1+x)/2, "single-step contrast",sign/color,"finite coefficient",s.str(12))
        else:
            assert to_arb(sign/color)<target
    assert passed==[Q(14,17),Q(16,17)]


if __name__=="__main__":
    main()
