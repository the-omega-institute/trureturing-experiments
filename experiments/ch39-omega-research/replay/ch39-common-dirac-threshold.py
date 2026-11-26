"""Equal one-loop offsets, with exact phase alignment and the same zero."""
from pathlib import Path
from runpy import run_path

from flint import arb, acb, ctx

ROOT = Path(__file__).parent
SELECT = run_path(str(ROOT/"ch39-rounded-phase-selection.py"))
MATCH,a = SELECT["MATCH"],SELECT["a"]


def masses(d,h,gamma,gain):
    ell2 = ((gamma/gain)**2+a("1/4")).sqrt()
    b = d*acb(((ell2+a("1/2"))/2).sqrt(), ((ell2-a("1/2"))/2).sqrt())
    c,e = b+d,b-d
    w,p = -h*c/d,-h*b/d
    assert (p-w-h).contains(0) and (p*c-w*b).contains(0)
    chi = c*e/b**2
    assert ((chi-1)/(chi+1)-acb(0,gain/(2*gamma))).contains(0)
    logs = [abs(m).log() for m in (p,w,c,b,e)]
    weights = [("1/5","3/10","4/5","1/10","3/5"),
               ("0","1/2","0","3/2","0"),("1/2","0","1/2","1","0")]
    offsets = [-a("4/3")*sum(a(t)*m for t,m in zip(row,logs)) for row in weights]
    scale = abs(w*b**3)
    assert all((v+a("2/3")*scale.log()).contains(0) for v in offsets)
    return offsets,scale


def verify():
    ctx.prec = 160
    ratio = MATCH["RATIO"]
    polynomials = ratio["product_polynomials"]()
    for n,kind,r,v,state in MATCH["matched_points"]():
        if kind!="center":
            continue
        z = 1/state[0]
        gamma = acb.zeta_zero(n).imag
        d,h = r*(a("13/10")*z).sqrt(),r*(a("7/2")*z).sqrt()
        epsilon = (state[5]/state[2]).sqrt()
        x = state[2].sqrt()
        offsets,scale = masses(d,h,gamma,624)
        change = (-scale.log()/9).exp()
        scores = MATCH["SCATTER"]["scalar_thresholds"](polynomials,r,a("1/12"))
        delta = MATCH["threshold_correction"](r*r,a("1/256"))
        low2,v2,lx2 = ratio["MODEL"]["matched"]([s+t+c for s,t,c in zip(scores,delta,offsets)],
                                                  epsilon,a("8/5"),r)
        assert (v2-v).contains(0) and (low2-1/(state[0]+3*state[1]).sqrt()).contains(0)
        assert (lx2.exp()-x*change).contains(0) and abs(change-1)<a(".00015")
        # All finite K>=1: t/d=L decreases to 1/sqrt(2) and |C|/d increases with L.
        Lmin = a("1/2").sqrt()
        Lmax = (gamma**2+a("1/4")).sqrt().sqrt()
        def cmod(L):
            return (L*L+1+(2*L*L+1).sqrt()).sqrt()
        small = (h/d)*d**4*Lmin**3*cmod(Lmin)
        large = (h/d)*d**4*Lmax**3*cmod(Lmax)
        assert x*(-small.log()/9).exp() < (4*arb.pi()).sqrt()
        assert epsilon*x*(-large.log()/9).exp()*v>1000
        # |E|=d L^2/cmod(L), cmod(L)<=L+1. All the other masses are larger
        # than this same elementary lower bound or its h/d multiple.
        assert min(d.lower(),h.lower())*Lmin**2/(Lmin+1)*v > 1000
        print(n,"fractional X,Y change",(change-1).str(9),
              "common offset in 2*pi/alpha units",offsets[0].str(9))


if __name__ == "__main__":
    verify()
