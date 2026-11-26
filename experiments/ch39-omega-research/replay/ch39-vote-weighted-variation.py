"""Fixed 49-height fiber for the vote-weighted mean, retaining the exact moment cancellations."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

root=Path(__file__).resolve().parent
tools=run_path(str(root/"ch39-fixed-poles-exposure.py"))
base,B,ball=tools["a"],tools["B"],tools["ball"]


def fiber():
    ctx.prec=480
    records=base["data"]("ch39-cutoff-threshold-inverses.json")
    original,changed=[base["exact_fit"](records[i]) for i in (2,1)]
    ids=sorted(set(original)|set(changed))
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    poles=[]
    for i in ids:
        h=prefix[i-1][0]
        poles.append((arb(h).union(arb(h+1))/arb(2)**400)**2)
    knots=[-j*j for j in base["OFFSETS"]]
    centers=[45,81,127,210,473]
    nodes=[poles[ids.index(i)] for i in centers]
    coeffs=[]
    for i,c in zip(centers,nodes):
        term=ball(original[i]-changed[i])*prod(1-c/v for k,v in zip(ids,poles) if k!=i)
        coeffs.append(term/(prod(c-k for k in knots)*prod(c-d for j,d in zip(centers,nodes) if j!=i)))
    def difference_without_y(y):
        poly=sum(v*prod(y-d for j,d in enumerate(nodes) if j!=i) for i,v in enumerate(coeffs))
        return poly*prod(y-j for j in knots[1:])/prod(1-y/c for c in poles)
    weights=[[ball(state.get(i,B(0))) for i in ids] for state in (changed,original)]
    ends=[[(1-t)*l+t*r for l,r in zip(*weights)] for t in (arb(1)/4,arb(3)/4)]
    return poles,ends,difference_without_y


def main():
    poles,ends,difference_without_y=fiber()
    muprime=difference_without_y(arb(0))
    p,modes,q02,_=base["physical_inputs"]()
    p,q02=ball(p),ball(q02)
    modes=list(map(ball,modes))
    # Linear expressions are enclosed by evaluating both endpoints first.
    mus=[sum(w/c for w,c in zip(row,poles)) for row in ends]
    mu=mus[0].union(mus[1])
    mass=arb(50772311235760).union(arb(50772311235762))
    q=q02/mass**2
    Z=(125/mass)**2
    count=sum(mu*c<q for c in poles)
    assert count==2 and mu*poles[count]>q
    Dt=-11*25*count*muprime/mu
    Nm=arb(25*count)
    for lam,n in zip(modes,base["MODES"]):
        def f(y):
            vals=[u*y-Z-lam+lam*(1-p)*sum(w/(1-y/c) for w,c in zip(row,poles)) for u,row in zip(mus,ends)]
            return vals[0].union(vals[1])
        for j in range(count+1):
            start=poles[j-1].upper() if j else arb(0)
            stop=poles[j].lower()
            l,r=start,stop
            for _ in range(82):
                m=(l+r)/2
                if f(m)<0:l=m
                else:r=m
            lo=l
            l,r=start,stop
            for _ in range(82):
                m=(l+r)/2
                if f(m)>0:r=m
                else:l=m
            hi=r
            assert f(lo)<0 and f(hi)>0
            y=lo.union(hi)
            assert start<y and y<stop and mu*y<q
            slopes=[u+lam*(1-p)*sum(w*c/(c-y)**2 for w,c in zip(row,poles)) for u,row in zip(mus,ends)]
            fy=slopes[0].union(slopes[1])
            Dt+=11*n*((muprime+lam*(1-p)*difference_without_y(y))/fy-muprime/mu)
            Nm+=n*(1-Z/(y*fy))
    tangent=mass*Dt/(22*Nm)
    assert arb("1.07130e-29")<tangent and tangent<arb("1.07132e-29")
    print("Weighted mean change",muprime.str(14),"parent mass derivative",tangent.str(16),"GeV.")
    return tangent


if __name__=="__main__":
    main()
