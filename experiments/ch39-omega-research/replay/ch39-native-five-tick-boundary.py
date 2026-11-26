"""First incompatible consecutive target for the fixed five-tick prime clock."""
from math import prod
from pathlib import Path
from runpy import run_path
import json

from flint import acb, arb, arb_mat, arb_poly, ctx


def separate(cert,ticks,precision,target_bound,pointwise,torus_bound,minimum_time_bits):
    root=Path(__file__).resolve().parent

    prior=run_path(str(root/"ch39-rh-shell-first-conflict.py"))
    ctx.prec,ctx.cap=2048,3
    ids=cert["nodes"]
    deg=len(ids)-1
    assert deg>=1 and len(cert["values"])==len(ids)+4
    nodes=[1/acb.zeta_zero(i).imag**2 if i else arb(0) for i in ids]
    h=[arb(x)/2**precision for x in cert["values"][:-4]]
    c=[arb(x)/2**precision for x in cert["values"][-4:]]
    levels=prior["read"](deg+1)
    V=arb_mat([[1/(1+j*j*z) for j in range(deg+1)] for z in nodes])
    coeff=V.solve(arb_mat([[v] for v in h])).entries()
    target=sum(a*f for a,f in zip(coeff,levels))
    assert target>arb(target_bound)
    cost=sum(abs(v).upper() for v in coeff)
    assert sum(abs(v).upper() for v in c)<2
    # These tolerances are on every time and phase observation respectively.
    assert target-arb(pointwise)>cost*arb(2)**(-minimum_time_bits)+2*arb(2)**-16
    def denom(z):
        return prod((1+j*j*z for j in range(1,deg+1)),start=arb(1))
    factors=[value*denom(z)/prod(z-o for j,o in enumerate(nodes) if i!=j)
             for i,(z,value) in enumerate(zip(nodes,h))]
    by_id=dict(zip(ids,h))
    log=[ticks*arb(x).log() for x in (2,3,arb(3)/2,6)]
    with ctx.workprec(288):
        gs=[x.imag for x in acb.zeta_zeros(1,3200)]
        assert gs[-1]<arb(2)**15
        for i,g in enumerate(gs,1):
            if i in by_id:
                r=by_id[i]
            else:
                z=1/g**2
                diffs=[z-o for o in nodes]
                r=prod(diffs)*sum(d/v for d,v in zip(factors,diffs))/denom(z)
            # The native law differs by less than 2^-96 in each component
            # throughout this height range, including the dyadic floor.
            assert r+sum(a*(g*l).cos() for a,l in zip(c,log))+2*arb(2)**-96 < arb(pointwise)
        end=1/gs[-1]**2
    x=arb_poly([0,1]); X=2*x-1
    c1,c2,c3,c4=c
    if torus_bound is None:
        torus=sum(abs(v).upper() for v in c)
    else:
        torus=arb(torus_bound)
        assert torus>abs(c1)
        prior["bernstein_test"]((torus-c1*X)**2-(c2+(c3+c4)*X)**2-(c3-c4)**2*(1-X**2),2)
    denominator=prod((1+j*j*end*x for j in range(1,deg+1)),start=arb_poly([1]))
    numerator=sum((d*prod((end*x-o for k,o in enumerate(nodes) if k!=i),start=arb_poly([1]))
                   for i,d in enumerate(factors)),arb_poly([0]))
    prior["bernstein_test"]((arb(pointwise)-torus)*denominator-numerator,deg)
    print("Tick",ticks,"separated at offsets 0 through",deg,"target",target.str(6),"pointwise upper",pointwise)


def main():
    root=Path(__file__).resolve().parent
    cert=json.loads((root/"ch39-native-five-tick-boundary.json").read_text())
    assert len(cert["nodes"])==54 and len(cert["support"])==57
    run_path(str(root/"ch39-prime-clock-height-horizon.py"))["main"]()
    separate(cert,5,48,".00076389",".0001",".436809",272)
    run_path(str(root/"ch39-native-clock-held-out.py"))["main"]([(cert["support"],53)],arb(0))


if __name__=="__main__":
    main()
