"""Calibrate the six rational bins, then test the threshold against input rounding."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_series, ctx

record = run_path(str(Path(__file__).with_name("ch39-finite-rh-pivot-mass.py")))


def main():
    ctx.prec,ctx.cap = 256,3
    masses = []
    for j in range(4):
        s = arb_series([arb(1)/2+j,1])
        xi = s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
        masses.append(xi[2]/xi[0] if j==0 else xi[1]/(2*j*xi[0]))
    for j in range(1,4):
        for i in range(j):
            h = (i*masses[i]+j*masses[j])/((i+j)*(masses[i]*masses[j]).sqrt())
            assert (h*record["den"]).floor().unique_fmpz()==record["bins"][j][i]
    h = 125**2
    a,b = arb(8).log()/2,arb(24).log()/2
    eig = [arb(6)]*7+[arb(4)]*6
    for n,k in ((3,5),(2,3)):
        eig += [3-arb(k).sqrt()]*n+[3+arb(k).sqrt()]*n
    def running(ms,couplings):
        inverse = [8*arb.pi()**2/arb(v)**2 for v in couplings]
        K = 5*inverse[0]-3*inverse[1]-2*inverse[2]
        F = arb(109)/15*(inverse[1]-inverse[2])-arb(23)/6*(inverse[0]-inverse[1])
        Z = sum((val+arb(h)/ms**2).log() for val in eig)/2
        A = K-44*(arb(ms)/1000).log()+12*a+4*b
        B = 15*F+22*Z-402*a-46*b
        ell,w = (448*A+16*B)/27808,(506*A-44*B)/27808
        assert w>b and ell>2 and arb(ms)**2*min(eig)+h>10**6
        return (ell-w).exp(),ms*ell.exp()
    # A point strictly inside the half-last-decimal input box.
    shifted = (".46745251",".6381149","1.058251")
    assert running(1137,shifted)[0]<arb("265.97205") and running(1138,shifted)[0]>arb("265.97207")
    central = (".467453",".63811","1.0583")
    lo,hi = arb("2034928.14"),arb("2034928.16")
    assert running(lo,central)[0]<arb("528.89869680") and running(hi,central)[0]>arb("528.89869682")
    _,crossing = running(lo.union(hi),central)
    assert crossing>arb("3.8981051e12") and crossing<arb("3.8981054e12")
    print("All six RH Gram entries certified inside the stored 96-bit bins.")
    print("Rounded-input three-node crossing bracket in scalar GeV: 1137 to 1138.")
    print("Minimum-residual middle node: scalar scale 2034928.14 to 2034928.16 GeV at the central gauge inputs.")


if __name__=="__main__":
    main()
