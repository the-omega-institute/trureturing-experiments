"""Minimum binomial degree of a shared phase law on the four truncated shell graphs."""
from itertools import combinations
from math import comb
from pathlib import Path
from runpy import run_path

home=Path(__file__).resolve().parent
original=run_path(str(home/"ch39-shell-local-nine-flux.py"))
model=run_path(str(home/"ch39-shell-third-shell.py"))


def solve(rows, values, modulus):
    """Solve a rectangular system over Z/(3^k), keeping a witness when solvable."""
    n=len(rows[0])
    a=[[v%modulus for v in row]+[val%modulus] for row,val in zip(rows,values)]
    perm=list(range(n))
    stages=[]
    start=0
    mod=modulus
    while mod>1 and start<len(a):
        pivot=next(((i,j) for i in range(start,len(a)) for j in range(start,n)
                    if a[i][j]%3),None)
        if pivot is None:
            if any(row[-1]%3 for row in a[start:]):
                return None
            for i in range(start,len(a)):
                a[i]=[v//3 for v in a[i]]
            mod//=3
            continue
        i,j=pivot
        a[start],a[i]=a[i],a[start]
        for row in a:
            row[start],row[j]=row[j],row[start]
        perm[start],perm[j]=perm[j],perm[start]
        inv=pow(a[start][start],-1,mod)
        a[start]=[inv*v%mod for v in a[start]]
        for i in range(start+1,len(a)):
            b=a[i][start]
            a[i]=[(x-b*y)%mod for x,y in zip(a[i],a[start])]
        stages.append(mod)
        start+=1
    y=[0]*n
    for i in reversed(range(start)):
        y[i]=(a[i][-1]-sum(v*w for v,w in zip(a[i][i+1:n],y[i+1:])))%stages[i]
    witness=[0]*n
    for i,j in enumerate(perm):
        witness[j]=y[i]
    assert all((sum(v*w for v,w in zip(row,witness))-rhs)%modulus==0 for row,rhs in zip(rows,values))
    return witness


def system(n,mod,sign,a):
    unit=mod//9
    rows=original["build"](n)
    labels=((2,1),(4*sign,6*sign))
    vals=[unit*labels[i][j] if j>=0 else 0 for _,_,i,j in original["pairs"]]
    for r,vs in model["fattest"](6).items():
        d0,d1,d2=model["shells"][r][0]
        ds={r-model["dot"](v,w) for v,w in combinations(vs,2)}
        assert sorted(d for d in ds if d<=d2)==[1,3,4,5,d2]
        for d in sorted(ds):
            if d==d0 or d>d2:
                continue
            rows.append([comb(d,k)-comb(d0,k) for k in range(1,n+1)])
            vals.append(a if d==d1 else unit*(1 if r==26 else 2) if d==d2 else 0)
    return rows,vals


def minimum(mod,sign,a):
    lo,hi=0,128
    assert solve(*system(hi,mod,sign,a),mod) is not None
    while hi-lo>1:
        mid=(hi+lo)//2
        if solve(*system(mid,mod,sign,a),mod) is not None:
            hi=mid
        else:
            lo=mid
    return hi


def main():
    for sign in (1,-1):
        degrees=[minimum(9,sign,a) for a in range(9)]
        order27=minimum(27,sign,19)
        assert degrees==[82,82,54,82,82,54,82,82,54] and order27==82
        print("nineteen-site orientation",sign,"degrees for nine phases",degrees,"robust order 27",order27)
    # Lucas at degrees <=81 forces F(5)=F(86) mod 3.  The degree-81
    # coefficient vanishes mod 3 by the unmarked pair 34,196 separated by 2*81.
    ds=original["sets"][1]
    assert {6,34,86,196}<=set(ds)
    for k in range(1,82):
        assert (comb(86,k)-comb(5,k) - (k==81))%3==0
        assert (comb(196,k)-comb(34,k) - 2*(k==81))%3==0


if __name__=="__main__":
    main()
