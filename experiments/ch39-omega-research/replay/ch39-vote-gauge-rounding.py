"""Gauge input rounding in the 119-mode chamber of the vote-weighted physical model."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

root=Path(__file__).resolve().parent
base=run_path(str(root/"ch39-cutoff-thresholds.py"))
ball=run_path(str(root/"ch39-fixed-poles-exposure.py"))["ball"]
N=119


def gauge(gs):
    d=[8*arb.pi()**2/g**2 for g in gs]
    return (sum(v*x for v,x in zip((355,-345,-10),d)),
            1000*((5*d[0]-3*d[1]-2*d[2])/44).exp())


def main():
    ctx.prec=320
    central=list(map(arb,(".467453",".63811","1.0583")))
    half=list(map(arb,(".0000005",".000005",".00005")))
    lower=[g+s*h for g,h,s in zip(central,half,(1,-1,-1))]
    upper=[g+s*h for g,h,s in zip(central,half,(-1,1,1))]
    V,_=gauge(central)
    vlo,qlo=gauge(lower);vhi,qhi=gauge(upper)
    # Inputs appear in both Q and V with the same monotonicity directions.
    assert qlo>arb("4.9987e13") and qhi<arb("5.0049e13") and vlo<V<vhi
    samples=base["data"]("ch39-vote-weighted-roots.json")
    a,b=[r["mass"] for r in samples]
    assert (a,b)==(50772311235760,50772311235762)
    lo,hi=50745000000000,50800000000000
    # Use the full enclosure of the certified reference weights to bound mu.
    rec=base["data"]("ch39-cutoff-threshold-inverses.json")[1]
    w=base["exact_fit"](rec)
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    def g2(i):
        h=prefix[i-1][0]
        return (arb(h).union(arb(h+1))/arb(2)**400)**2
    # The interpolation proof bounds the deviation of every admissible mean.
    mu=sum(ball(u)/g2(i) for i,u in w.items())+arb(0,arb("4e-55"))
    open_y=[g2(1),g2(2)]
    for block in samples[0]["roots"]:
        open_y.extend(arb(l).union(arb(h))/arb(2)**64 for l,h in block)
    # For H(M)=M^2 H0 + Hh, with Hh positive, the eigenvalue at M is
    # between min(1,(M/a)^2) and max(1,(M/a)^2) times the one at a.
    assert all(lo**2*mu*y>10**19 and hi**2*mu*y<qlo**2 for y in open_y)
    assert lo**2*mu*g2(3)>qhi**2 and qhi**2<8*lo**2
    # The remaining mixed and pure modes are above the third pole.
    # Hellmann-Feynman, or min-max difference quotients, gives a slope in
    # [2(1-125^2/10^19),2] for each open log squared mass vs log M.
    correction=(1-arb(125**2)/10**19).union(arb(1))
    for corner,lv,hv in [(vlo,50745372380748,50745372380753),(vhi,50799263747305,50799263747311)]:
        result=arb(a).union(arb(b))*((corner-V)/(22*N*correction)).exp()
        assert lo<lv<result and result<hv<hi
        print("Mass at gauge rounding corner in GeV:",lv,hv)
    print("All gauge inputs in the half-last-digit box stay in the 119-mode chamber.")


if __name__=="__main__":
    main()
