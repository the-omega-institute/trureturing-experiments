"""Classify all exact zero fields of the quaternion four-cycle heap comparison."""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

directory=Path(__file__).parent
heap=run_path(str(directory/"ch39-rh-heap-finite-shell-metric.py"))
mul,conj=heap["mul"],heap["conj"]


def compose(p,q):
    return tuple(p[i] for i in q)


def parity(p):
    return (-1)**sum(a>b for a,b in combinations(p,2))


def quotient():
    states=list(permutations(range(4)))
    klein=[p for p in states if p==states[0] or (parity(p)==1 and all(i!=v for i,v in enumerate(p)))]
    four=[p for p in states if parity(p)<0 and all(i!=v for i,v in enumerate(p))]
    cosets=sorted({tuple(sorted(compose(p,h) for h in klein)) for p in states})
    cid={p:i for i,c in enumerate(cosets) for p in c}
    even=[i for i,c in enumerate(cosets) if parity(c[0])==1]
    odd=[i for i,c in enumerate(cosets) if parity(c[0])==-1]
    return states,klein,four,cosets,cid,even,odd


def main():
    states,klein,four,cosets,cid,even,odd=quotient()
    edges=list(product(even,odd))
    patterns={}
    for mask in range(1<<len(edges)):
        orth={edge for i,edge in enumerate(edges) if mask&(1<<i)}
        if any(sum(v in e for e in orth)%2 for v in range(6)):
            continue
        components=list(range(6))
        for left,right in set(edges)-orth:
            old,new=components[right],components[left]
            components=[new if c==old else c for c in components]
        if any(components[left]==components[right] for left,right in orth):
            continue
        labels={v:i for i,v in enumerate(sorted(set(components)))}
        patterns[mask]=[labels[v] for v in components]
    assert len(patterns)==7
    basis=((1,0,0,0),(0,1,0,0),(0,0,1,0))
    trails=[(p,compose(p,g),compose(p,compose(g,g))) for p,g in product(states,four)]
    quaternions=[(tuple(s*v for v in b)) for b,s in product(((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)),(-1,1))]
    predicted={(v,w):mul(mul(w,conj(v)),w) for v,w in product(quaternions,repeat=2)}
    frame_classes=0
    scalar_layers={edge:0 for edge in trails}
    for mask,labels in patterns.items():
        nlines=len(set(labels))
        assert (mask==0 and nlines==1) or (mask.bit_count()==6 and nlines==3)
        signs={}
        for i,c in enumerate(cosets):
            p=c[0]
            signs[p]=1
            for g in four:
                q=compose(p,g)
                r=compose(q,g)
                value=1 if labels[i]==labels[cid[q]] else -1
                assert r not in signs or signs[r]==value
                signs[r]=value
            assert all(h in signs for h in c)
        for flips in product((-1,1),repeat=6):
            field={p:tuple(flips[cid[p]]*signs[p]*v for v in basis[labels[cid[p]]]) for p in states}
            assert all(predicted[field[p],field[q]]==field[r] for p,q,r in trails)
        frame_classes+=2**(6-nlines)
        if mask:
            for edge in trails:
                p,q,r=edge
                # Sum the number of body axes fixed by conjugation, using three
                # cyclic axes as a discrete isotropic ensemble.
                scalar_layers[edge]+=3 if labels[cid[p]]==labels[cid[q]] else 1
    assert frame_classes==80 and set(scalar_layers.values())=={10}
    take=Q(10,18)
    assert 16*take==Q(80,9) and (16*take)/Q(8,3)==Q(10,3)
    ctx.prec,ctx.cap=1536,3
    arc=heap["arc"]
    data=arc["fixed"]["inputs"]()
    first=arc["argument_for"](arb(33)/34,data)
    d=arc["argument_for"](arb(14)/17,data)/(3*first)
    stiffness=(d-1)/(arb(80)/9-arb(8)/3*d)
    assert 3<d<arb(10)/3 and stiffness>0
    normalizer=1+arb(8)/3*stiffness
    print("Exact projective patterns 1+6; labelled configurations up to a common frame: 32+48=80.")
    print("Nonconstant-pattern Hessian average: 5/9 three-point and 4/9 two-point terms.")
    print("Stiffness",stiffness.str(12),"neutral response",arc["response"](3*first/normalizer,data).str(13))


if __name__=="__main__":
    main()
