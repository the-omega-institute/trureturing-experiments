"""A finite neighborhood of unit kinetic scale preserves the exact five-tick system."""
from pathlib import Path
from runpy import run_path

from flint import acb,arb,arb_mat,ctx

prior=run_path(str(Path(__file__).with_name("ch39-rh-shell-smallest-tick.py")))


def main():
    ctx.prec,ctx.cap=1536,3
    heights=[acb.zeta_zero(i).imag for i in prior["support"]]
    poles=[1/g**2 for g in heights]
    logs=[arb(v).log() for v in (2,3,arb(3)/2,6)]
    matrix=arb_mat([[1/(1+j*j*z) for z in poles] for j in range(42)]+
                   [[(5*g*log).cos() for g in heights] for log in logs])
    target=prior["prior"]["read"](42)+[arb(0)]*4
    inv=matrix.inv()
    weights=(inv*arb_mat([[v] for v in target])).entries()
    assert min(weights)>arb(".00079")
    # Delta M = (scale - 1) times the displayed divided difference. Form the
    # preconditioned matrix before taking a norm, preserving cancellations.
    for bits in (180,160,140,120,100,80):
        rad=arb(2)**(-bits)
        factor=1+arb(0,rad)
        divided=arb_mat([[-j*j*z/((1+j*j*z)*(1+factor*j*j*z)) for z in poles]
                         for j in range(42)]+[[arb(0)]*46 for _ in logs])
        E=inv*divided
        norm=rad*max(sum(abs(E[i,j]) for j in range(46)) for i in range(46))
        if norm<1 and norm*max(weights)/(1-norm)<min(weights):
            best=bits
    assert best<=180
    print("All 42 exact readings and four five-tick angle constraints have positive weights for |scale - 1| <= 2^-",best)


if __name__=="__main__":
    main()
