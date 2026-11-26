"""Two-parameter determinant of the kinetic rows on the two-dimensional moment kernel."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
tom=run_path(str(root/"ch39-two-probe-tomography.py"))
s=tom["s"]
base,B,rounded,multiply=s["base"],s["B"],s["rounded"],s["multiply"]
unit=10**6


def numerator(tri,k):
    a,b=tri.fits[k],tri.fits[2]
    ids=sorted(set(a)|set(b))
    assert len(ids)==49
    missing,=set(tri.ids)-set(ids)
    r=tri.poles
    nodes=(45,81,127,210,473)
    terms=[]
    for i in nodes:
        v=r[i]
        top=multiply(a[i]-b[i],*(r[j]-v for j in ids if j!=i))
        den=multiply(*(B(-j*j)-v for j in base["OFFSETS"] if j),*((r[j]-v)/B(unit) for j in nodes if j!=i))
        row=[rounded(top/den)]
        for j in nodes:
            if j!=i:row=tom["times_linear"](row,r[j]/B(unit))
        terms.append(row)
    poly=[rounded(sum(row,B(0))) for row in zip(*terms)]
    return tom["times_linear"](poly,r[missing]/B(unit))


def main():
    tri=s["ThreeFits"]()
    a,b=numerator(tri,0),numerator(tri,1)
    # (t^i s^j - t^j s^i)/(t-s) for i>j is a finite geometric sum.
    bez=[[B(0) for _ in range(5)] for _ in range(5)]
    for i in range(1,6):
        for j in range(i):
            v=rounded(a[i]*b[j]-a[j]*b[i])
            for h in range(i-j):
                bez[i-1-h][j+h]=rounded(bez[i-1-h][j+h]+v)
    assert all(v.lo>0 for v in a)
    assert all(v.lo>0 for v in b[:4]) and all(v.hi<0 for v in b[4:])
    assert (a[0]/b[0]).hi<(a[1]/b[1]).lo
    q=Q(1677,1000)
    def value(poly,x):
        out=B(0)
        for v in reversed(poly):out=rounded(out*B(x)+v)
        return out
    level=[v-B(q)*w for v,w in zip(a,b)]
    points=[Q(v,1000) for v in (5,15,32,55)]
    for u,sign in zip(points,(-1,1,-1,1)):
        assert (B(sign)*value(level,u)).lo>0
    # Put t=1/r, and bound B(t*v,v) on strips in 0<=t<=1/5.
    for m in range(96):
        lo,hi=Q(m,480),Q(m+1,480)
        poly=[B(0) for _ in range(9)]
        for i in range(5):
            for j in range(5):
                poly[i+j]+=bez[i][j]*B(lo**i,hi**i)
        positive_orthant(poly)
    for multiple,interior in [(2,Q(".015")),(3,Q(".012")),(4,Q(".0095"))]:
        def scaled_det(x):
            return sum((bez[i][j]*B(x**(i+j)*multiple**j) for i in range(5) for j in range(5)),B(0))
        for u,sgn in ((Q(".004"),1),(interior,-1),(Q(".06"),1)):
            assert (B(sgn)*scaled_det(u)).lo>0
    print("Four arbitrary distinct fresh real spacings suffice; three can be redundant.")
    print("Any pair of fresh spacings with ratio at least five suffices for all positive known scales.")
    print("Ratios 2, 3 and 4 each have exceptional positive settings.")


def positive_orthant(poly):
    n=len(poly)-1
    initial=[v.lo/comb(n,j) for j,v in enumerate(poly)]
    stack=[(initial,0,True)]
    while stack:
        c,depth,right_unbounded=stack.pop()
        if min(c)>=0 and c[0]>0 and (c[-1]>0 or right_unbounded):continue
        assert depth<25
        left,right=[c[0]],[c[-1]]
        while len(c)>1:
            c=[(x+y)/2 for x,y in zip(c,c[1:])]
            left.append(c[0])
            right.append(c[-1])
        stack.extend([(left,depth+1,False),(list(reversed(right)),depth+1,right_unbounded)])


if __name__=="__main__":
    main()
