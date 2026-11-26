"""Sharp time-row count for identifying the 50-site support on the ordered RH grid."""
from fractions import Fraction as Q
from itertools import groupby
from math import prod
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
lib=run_path(str(root/"ch39-vote-weighted-integers.py"))
base,B,rounded,multiply,absolute=[lib[x] for x in ("base","B","round_out","multiply","absolute")]
spread=(1445,1477,1510,1539,1687,1725,2058,2249,2511,2805,3274,3739,4364,4874,
        5566,6497,6642,7256,7754,8471,9461,11044,13180,16441,21436,32624,73918)


def roster():
    records=base["data"]("ch39-cutoff-threshold-inverses.json")
    ids=sorted(set.union(*(set(r["support"]) for r in records)))
    runs=[[v for _,v in group] for _,group in groupby(enumerate(ids),lambda p:p[1]-p[0])]
    roots=ids+[v[-1] for v in runs[1:] if len(v)%2]
    # One barycentric signed dependence suffices for any <=70 distinct rows.
    path=[]
    for i in range(1,ids[-1]+2):
        if len(path)%2==0 or i in ids:path.append(i)
    assert len(ids)==50 and len(roots)==70 and len(path)==71
    assert set(path[1::2])<=set(ids) and set(path[::2])-set(ids)
    # Between neighboring roots, log(P/W) is strictly concave in log z.
    # d/dz [z d(log(P/W))/dz] is a sum of strictly negative terms.
    ends=set()
    for run in runs:
        if run[0]>1:ends.add(run[0]-1)
        ends.add(run[-1]+1)
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    z={i:rounded(B(Q(2**800,(prefix[i-1][0]+1)**2),Q(2**800,prefix[i-1][0]**2))) for i in set(ids)|ends}
    return roots,ends,z


def condition(extra,exponent,worst):
    roots,ends,z=roster()
    A=sorted([j*j for j in base["OFFSETS"] if j]+[j*j for j in extra])
    assert len(A)==len(set(A))==70
    roots=[z[i] for i in roots]
    # Exact normalization cancels the coefficient of the constant row.
    cost=B(0)
    for a in A:
        numer=multiply(*(v+B(Q(1,a)) for v in roots))
        den=prod(abs(1-Q(k,a)) for k in A if k!=a)
        cost+=rounded(numer/B(den))
    gap={i:rounded(multiply(*(z[i]-v for v in roots))/multiply(*(B(1)+B(a)*z[i] for a in A))) for i in ends}
    assert all(v.lo>0 for v in gap.values())
    assert all(gap[worst].hi<v.lo for i,v in gap.items() if i!=worst)
    cond=cost/gap[worst]
    assert 2**(exponent-1)<cond.lo<=cond.hi<2**exponent
    print("70 nonconstant time rows, minimum at zero",worst,"leakage factor below 2^",exponent)


def main():
    condition(range(43,70),420,2322)
    condition(spread,245,85)
    print("71 time rows including normalization suffice, and any 70 admit a positive off-roster alternative.")


if __name__=="__main__":
    main()
