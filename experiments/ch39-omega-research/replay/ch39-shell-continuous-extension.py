"""A continuous first reference phase on the two tied twelve-site shells."""
from functools import lru_cache
from pathlib import Path
from runpy import run_path
from flint import acb, arb

model = run_path(str(Path(__file__).with_name("ch39-shell-third-shell.py")))
penalty, choices = model["penalty"], model["choices"]
measure, distance = model["measured"], model["distance"]
pi = arb.pi()


def point(h):
    base = sum(distance(measure(h, r), ref) for r,ref in model["original"].items())
    return base + penalty*h*h, {r: measure(h, r) for r in model["extended"]}


@lru_cache(None)
def zero(n):
    return point(acb.zeta_zero(n).imag)


def refs(x):
    return {r: model["readout"](model["shells"][r][1], pi*x/9, pi*j/9)
            for r,j in ((26,1),(29,2))}


def score(layers, state, ref):
    base, obs = state
    return base + sum(distance(obs[r], ref[r]) for r in layers)


def delta(layers, n, x):
    # The reference norm drops out before interval evaluation.
    ref = refs(x)
    old, ys = zero(719)
    other, zs = zero(n)
    return other-old + sum((z-y)*(z+y-2*t) for r in layers
                           for z,y,t in zip(zs[r],ys[r],ref[r]))


def coefficients(layers, n):
    # x is in units pi/9. Each squared deficit difference has only one harmonic.
    f0, f90, f180 = (delta(layers,n,x) for x in (0,arb(9)/4,arb(9)/2))
    constant = (f0+f180)/2
    return constant, f0-constant, f90-constant


def harmonic(coef, x):
    s,c = (2*pi*x/9).sin_cos()
    return coef[0]+coef[1]*c+coef[2]*s


left,right = arb(31)/5,arb(639)/100
slab = left.union(right)


def contenders():
    # Cover every zero whose energy interval could meet the target in this slab.
    ref = refs(slab)
    target = {rs:score(rs,zero(719),ref) for rs in choices}
    cap = max(int((s/penalty).sqrt().upper().ceil().unique_fmpz()) for s in target.values())
    assert all(penalty*cap**2 > s for s in target.values())
    found = {rs: set() for rs in choices}
    def scan(lo,hi,parts):
        state = point(lo.union(hi))
        keep = [rs for rs in parts if not score(rs,state,ref)>target[rs]]
        if not keep:
            return
        if hi-lo>arb(1)/4:
            mid=(hi+lo)/2
            scan(lo,mid,keep)
            scan(mid,hi,keep)
            return
        a,b=(v.zeta_nzeros().unique_fmpz() for v in (lo,hi))
        assert a is not None and b is not None
        for n in range(int(a)+1,int(b)+1):
            if n==719:
                continue
            for rs in keep:
                if not delta(rs,n,slab)>0:
                    found[rs].add(n)
    for h in range(cap):
        scan(arb(h), arb(h+1), choices)
    print("initial slab", left.str(8),right.str(8),"high tail",cap)
    return found


def root(co, increasing):
    s,c = (2*pi*slab/9).sin_cos()
    slope = -co[1]*s+co[2]*c
    sign = 1 if increasing else -1
    assert sign*slope > 0
    lo,hi = left,right
    assert sign*harmonic(co,lo)<0<sign*harmonic(co,hi)
    for _ in range(66):
        mid=(lo+hi)/2
        value=sign*harmonic(co,mid)
        if value<0:
            lo=mid
        else:
            assert value>0
            hi=mid
    return lo.union(hi)


def main():
    # These finitely many competitors suffice to exclude the entire rest of the circle.
    ids=(2447,2454,2619,357,651,858,583)
    witnesses=[coefficients(rs,n) for rs in choices[:2] for n in ids]
    def exclude(lo,hi):
        x=lo.union(hi)
        if any(harmonic(co,x)<0 for co in witnesses):
            return
        assert hi-lo>arb(1)/100000
        mid=(lo+hi)/2
        exclude(lo,mid)
        exclude(mid,hi)
    exclude(arb(0),left)
    exclude(right,arb(9))
    found=contenders()
    boundaries={((26,),651):True, ((29,),858):False}
    endpoints=[]
    for rs,ns in found.items():
        for n in sorted(ns):
            co=coefficients(rs,n)
            if (rs,n) in boundaries:
                endpoints.append(root(co,boundaries[rs,n]))
            else:
                # Harmonic representation avoids dependency in the individual squared deficits.
                assert harmonic(co,slab)>0
    assert all(n in found[rs] for rs,n in boundaries)
    lo,hi=endpoints
    assert lo<hi
    print("complete common interval, in units pi/9",lo.str(16),hi.str(16))
    print("additional phase in units pi", (lo/9).str(16),(hi/9).str(16))
    for m in (9,18,27):
        hits=[k for k in range(m) if lo<arb(9)*k/m<hi]
        assert hits==([19] if m==27 else [])
        print("common order",m,"angle numerators",hits)


if __name__=="__main__":
    main()
