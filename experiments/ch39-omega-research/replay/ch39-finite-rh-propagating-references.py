"""Propagating RH reference fields on the 48 shell edges and their one-loop crossing."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

record = run_path(str(Path(__file__).with_name("ch39-finite-rh-pivot-mass.py")))


def gram(i,j):
    lo,hi = record["read"](i,j)
    to_ball = lambda q:arb(q.numerator)/q.denominator
    return to_ball(lo).union(to_ball(hi))


def root_bracket(poly,lo,hi):
    lo,hi = Fraction(lo),Fraction(hi)
    to_ball = lambda q:arb(q.numerator)/q.denominator
    increasing = poly(to_ball(lo))<0
    assert (poly(to_ball(hi))>0) == increasing and not poly(to_ball(hi)).contains(0)
    for _ in range(46):
        mid = (lo+hi)/2
        f = poly(to_ball(mid))
        if f.contains(0):
            break
        if (f<0)==increasing:
            lo = mid
        else:
            hi = mid
    return to_ball(lo).union(to_ball(hi))


def main():
    ctx.prec=224
    d = [8*arb.pi()**2/arb(v)**2 for v in (".467453",".63811","1.0583")]
    k = 5*d[0]-3*d[1]-2*d[2]
    f = arb(109)/15*(d[1]-d[2])-arb(23)/6*(d[0]-d[1])
    q = 1000*(k/44).exp()
    a,b,c = gram(1,0),gram(1,2),gram(0,2)
    det = 1+2*a*b*c-a*a-b*b-c*c
    pivot = det/(1-c*c)
    bands = [(arb(6),7),(arb(4),6)]
    for j,n in ((5,3),(3,2)):
        bands.extend([(3-arb(j).sqrt(),n),(3+arb(j).sqrt(),n)])
    def exposure(mf):
        y = q/mf
        spectra = [(1+c,25),(1-c,25)]
        for l,n in bands:
            def poly(t):
                return (t-l-(125/mf)**2)*((t-1)**2-c*c)-l*((a*a+b*b)*(t-1)+2*a*b*c)
            for left,right in [("0",".0001"),(".001",".01"),("1","9")]:
                spectra.append((root_bracket(poly,left,right),n))
        total = arb(0)
        active = 0
        for mass2,n in spectra:
            assert mass2>0 and mass2.sqrt()*mf>1000
            term = (y/mass2.sqrt()).log()
            if term>0:
                active+=n
                total+=n*term
            else:
                assert term<0
        assert sum(n for _,n in spectra)==119 and active==71
        assert y<arb(8).sqrt()
        return total
    lo,hi = arb("2.23915e14"),arb("2.23917e14")
    assert 22*exposure(lo)>15*f and 22*exposure(hi)<15*f
    # For y >= sqrt(8), cycle modes, interlacing and the Schur trial
    # vector bound the exposure without computing any large-scale roots.
    w0=arb(8).log()/2
    lower = 50*w0-arb(25)/2*(1-c*c).log()+23*(w0-arb(2).log()/2)
    lower -= arb(23)/2*(arb(1)/64+6*pivot/8).log()
    assert 22*lower>15*f and 22*73>402+46
    print("One-loop crossing Q",q.str(15),"GeV; 71 of 119 extra scalar doublets open, no vectorlike fermion open.")
    print("Propagating-reference parent mass between 2.23915e14 and 2.23917e14 GeV.")


if __name__=="__main__":
    main()
