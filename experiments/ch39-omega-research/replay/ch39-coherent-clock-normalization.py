"""Conditional kinetic normalization under exact and approximate shell-clock commutation."""
from fractions import Fraction as Q
from itertools import combinations
from math import isqrt
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
source = run_path(str(root/"ch39-phase-support-resolution.py"))
base, B, rounded = source["base"], source["B"], source["rounded"]


def root_cell(value):
    scale = 1 << 128
    assert value.lo > 0
    return B(Q(isqrt((value.lo*scale**2)//1),scale),
             Q(isqrt((value.hi*scale**2)//1)+1,scale))


def clock_pairs(ids):
    rows = base["observations"](ids)[0][-4:]
    den = base["SCALE"]**2
    for i,j in combinations(range(len(ids)),2):
        gaps = []
        for row in rows:
            l,h = row[i][0]-row[j][1], row[i][1]-row[j][0]
            gaps.append(B(Q(min(l*l,h*h) if l*h>0 else 0,den),Q(max(l*l,h*h),den)))
        yield i,j,gaps


def check(weights, expected, gram=None):
    ids = sorted(weights)
    n = len(ids)
    total = [B(0) for _ in expected]
    pairs = [[] for _ in expected]
    gain = energy = B(0)
    prefix = base["data"]("ch39-tick-separators.json")["prefix"]
    z = [source["amb"]["upper"]["height_cell"](prefix[i-1][0]) for i in ids]
    if gram:
        vectors = gram["rows"]
        den = 1 << (2*gram["bits"])
        assert len(vectors)==n and len({len(v) for v in vectors})==1
        assert all(Q(sum(k*k for k in v),den)<zi.lo for v,zi in zip(vectors,z))
    for i,j,terms in clock_pairs(ids):
        ww = weights[ids[i]]*weights[ids[j]]
        for k,depth in enumerate((1,2,4)):
            gap = rounded(sum(terms[:depth],B(0)))
            assert gap.lo > 0
            total[k] = rounded(total[k] + B(2)*ww/gap)
            pairs[k].append((gap,ids[i],ids[j]))
        if gram:
            a = B(Q(sum(k*l for k,l in zip(vectors[i],vectors[j])),den))
            gain = rounded(gain + B(2)*root_cell(ww)*a)
            energy = rounded(energy + B(2)*gap*a*a)
    for (low,high,i,j),v,gaps in zip(expected,total,pairs):
        assert Q(low)**2<v.lo<=v.hi<Q(high)**2
        g,u,vv=min(gaps,key=lambda row:row[0].lo)
        assert (u,vv)==(i,j) and all(g.hi<h.lo for h,ii,jj in gaps if (ii,jj)!=(i,j))
    if gram:
        m1=sum((weights[i]*x for i,x in zip(ids,z)),B(0))
        m2=sum((weights[i]*x**2 for i,x in zip(ids,z)),B(0))
        excess=m2/m1-m1
        assert excess.lo>0 and total[2].hi/2**18<excess.lo**2
        assert Q("0.00228")<gain.lo and excess.hi<gain.lo
        assert energy.hi<Q(1,2**16)
        print("Positive Gram witness at eight bits; nine bits exclude the second-moment response.")
    print(n,"sites; sharp coherence factors bounded above by",*[a[1] for a in expected])


def main():
    reference,moved,*_ = source["tables"]()
    gram = base["data"]("ch39-coherent-clock-gram.json")
    check(reference,[("63.7623","63.7624",4,11),("1.93306","1.93307",16,25),
                     (".93092",".93093",9,195)],gram)
    check(moved,[("81.9669","81.9670",39,453),("1.93203","1.93204",16,25),
                 (".93041",".93042",9,195)])


if __name__=="__main__":
    main()
