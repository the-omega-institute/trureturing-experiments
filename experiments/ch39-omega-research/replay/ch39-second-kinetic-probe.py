"""Optimal uniform error exponent for a second kinetic probe of the fixed-scale pair."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
s=run_path(str(root/"ch39-fixed-scale-common-reading.py"))
base,B=s["base"],s["B"]
GRID=2**182


def grid(cell):
    return cell.lo*GRID//1,-((-cell.hi*GRID)//1)


def product(a,b):
    corners=[x*y for x in a for y in b]
    return min(corners)//GRID,-((-max(corners))//GRID)


def magnitude(a):
    return max(abs(a[0]),abs(a[1]))


def local_bound(poles,weights,c,lo,hi):
    mid,rad=(lo+hi)/2,(hi-lo)/2
    totals=[[0,0] for _ in range(10)]
    remainder=0
    for v,w in zip(poles,weights):
        inv=B(1)/(v+B(mid))
        ratio=grid(B(rad)*inv)
        assert 0<=ratio[0]<=ratio[1]<GRID
        term=product(w,grid(inv))
        for coeff in totals:
            coeff[0]+=term[0]
            coeff[1]+=term[1]
            term=product(term,ratio)
        remainder-=(-magnitude(term)*GRID)//(GRID-ratio[1])
    # Powers of (mid-x)/rad range from -1 to 1.
    return Q(sum(map(magnitude,totals))+remainder,GRID)/c


def main():
    tri=s["ThreeFits"]()
    c=Q(base["data"]("ch39-fixed-scale-key.json")[0]["c"])
    _,models=tri.models(c)
    poles=[tri.poles[i] for i in tri.ids]
    delta=[r[1]-l[1] for l,r in zip(models[0].ref,models[1].ref)]
    eta=1024
    contrast=sum((w/(B(c)*v+B(eta)) for w,v in zip(delta,poles)),B(0))
    assert Q("-6.5581e-12")<contrast.lo<contrast.hi<Q("-6.5578e-12")<-Q(1,2**38)
    weights=[grid(w) for w in delta]
    cap=Q(1,2**37)
    tail=2**40
    assert Q(sum(map(magnitude,weights)),GRID)/(c*tail)<cap
    cuts=[0]+[2**k for k in range(7,41)]
    stack=[(Q(x),Q(y),0) for x,y in zip(cuts,cuts[1:])]
    passed=0
    while stack:
        lo,hi,depth=stack.pop()
        if local_bound(poles,weights,c,lo,hi)<cap:
            passed+=1
        else:
            assert depth<18
            mid=(lo+hi)/2
            stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
    print("At eta=1024, two disjoint bands of radius 2^-39.")
    print("Every nonnegative eta admits overlapping bands at radius 2^-38;",passed,"finite Taylor intervals and one tail.")


if __name__=="__main__":
    main()
