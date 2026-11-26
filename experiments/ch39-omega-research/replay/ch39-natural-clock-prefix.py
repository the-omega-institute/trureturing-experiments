"""How long the natural RH weights can be completed by unconstrained two-clock phases."""
import json
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
integers=run_path(str(root/"ch39-vote-weighted-integers.py"))
base,B,rounded=integers["base"],integers["B"],integers["round_out"]
cheb=run_path(str(root/"ch39-tick-separators-integers.py"))["chebyshev_upper"]


def calibrate():
    from flint import arb,arb_series,ctx
    ctx.prec,ctx.cap=640,3
    s=arb_series([arb(1)/2,1])
    xi=s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
    return int((xi[2]/xi[0]*2**288).floor().unique_fmpz())


def margins(total,observed):
    a,b,c,d=observed
    v=B(1)-total
    def square(x):
        lower=0 if x.lo<=0<=x.hi else min(x.lo*x.lo,x.hi*x.hi)
        return B(lower,max(x.lo*x.lo,x.hi*x.hi))
    return [v-c,v-d,v+c,v+d,(v-c)*(v-d)-square(a+b),(v+c)*(v+d)-square(a-b)]


def main():
    stored=base["data"]("ch39-natural-mass-cell.json")
    mass=B(Q(stored,2**288),Q(stored+1,2**288))
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    total=B(0)
    sums=[B(0)]*4
    seen=[None]*4
    assert len(prefix)==3200 and prefix[-1][0]+1<2**415
    for i,(h,*angles) in enumerate(prefix,1):
        w=rounded(B(Q(2**800,(h+1)**2),Q(2**800,h*h))/mass)
        total=rounded(total+w)
        for j,angle in enumerate(angles):
            l,r=cheb(angle,5)
            val=B(Q(l,2**72)-Q(1,2**96),Q(r,2**72)+Q(1,2**96))
            sums[j]=rounded(sums[j]+rounded(w*val))
            if seen[j] is None:
                norm=sums[j] if sums[j].lo>=0 else B(-sums[j].hi,-sums[j].lo) if sums[j].hi<=0 else B(0,max(-sums[j].lo,sums[j].hi))
                if norm.lo>1-total.lo:seen[j]=i
                else:assert norm.hi<1-total.hi
        if i in (43,44):
            tests=margins(total,sums)
            if i==43:assert min(t.lo for t in tests)>0
            else:assert min(t.hi for t in tests)<0
    assert seen==[None,83,270,85]
    wtail=B(1)-total
    complete=[v+B(-wtail.hi,wtail.hi) for v in sums]
    assert complete[1].hi<Q("-.134") and complete[2].hi<Q("-.060") and complete[3].lo>Q(".142")
    # Monotonicity under attenuation follows by adding a zero-mean torus measure.
    def feasible(factor):
        return margins(total*B(factor),[x*B(factor) for x in sums])
    low,high=Q(0),Q(1)
    for _ in range(40):
        mid=(low+high)/2
        t=feasible(mid)
        if min(v.lo for v in t)>0:low=mid
        else:
            assert min(v.hi for v in t)<0
            high=mid
    assert Q(".825047100802")<low<high<Q(".825047100804")
    print("First incompatible natural-weight prefix 44; separate first failures",seen)
    print("Maximal retained fraction of the first 3200 weights",float(low),float(high))
    print("Ordinary TV between the natural measure and any four-null clock measure exceeds 0.071.")


if __name__=="__main__":
    import sys
    path=root/"ch39-natural-mass-cell.json"
    if "--write" in sys.argv:path.write_text(json.dumps(calibrate())+"\n")
    if "--origin" in sys.argv:assert json.loads(path.read_text())==calibrate()
    run_path(str(root/"ch39-prime-clock-height-horizon.py"))["main"](5,15)
    main()
