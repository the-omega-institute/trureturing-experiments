"""Integer gains with aligned determinant phases over the full five-digit mass window."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb,arb,ctx

ROOT = Path(__file__).parent
MATCH = run_path(str(ROOT/"ch39-quartic-soft-adjoint-match.py"))
PLANCK, RATIO, a = MATCH["PLANCK"],MATCH["RATIO"],MATCH["a"]


def target_gap(gamma, gain, sign):
    ell = ((gamma/gain)**2 + a("1/4")).sqrt().sqrt()
    u = 1/ell
    stretch = (1+u*u+u*(2+u*u).sqrt()).sqrt()
    return a("13/35").sqrt().sqrt().sqrt()/(ell*(stretch**sign).sqrt().sqrt())


def matches():
    products = RATIO["product_polynomials"]()
    sigma,lam,k = a("1/256"),a("1/12"),a("8/5")
    qbox = a(".16").union(a(".17"))**2
    assert (a("127/28")/qbox+a("109/84")*sigma/(1+sigma*qbox)
            -a("97/7")/(2+qbox)+a("1/21")/(1+2*qbox))>0
    assert (-a("109/14")/qbox-a("107/42")*sigma/(1+sigma*qbox)
            +a("369/14")/(2+qbox)-a("1/21")/(1+2*qbox))<0
    def normalized(q):
        return q,(1+sigma*q)/(1+sigma),(q+2)/3,(1+2*q)/3
    def shift(q, powers):
        return sum(a(e)*f.log() for e,f in zip(powers,normalized(q))).exp()
    for n,gamma,eps in RATIO["high_zero_readings"]():
        base = MATCH["SCATTER"]["scalar_thresholds"](products,arb(1),lam)
        corr = MATCH["threshold_correction"](arb(1),sigma)
        l0,v0,_ = RATIO["MODEL"]["matched"]([s+ds for s,ds in zip(base,corr)],eps,k,arb(1))
        scale = (2*arb.pi()*(gamma**2+a("1/4")+30+a("2/15"))/(3*arb(5).sqrt()/2)).sqrt()
        def mass(r):
            return scale*v0/l0*shift(r*r,("-109/14","-107/42","369/14","-1/42"))
        def gap(r):
            return (a("13/10")*(1+3/k**2)).sqrt()*l0*r*shift(r*r,("113/28","109/84","-97/7","1/42"))
        def root(fn,target,increasing):
            return PLANCK["enclose_root"](fn,target,F(".16"),F(".17"),increasing)
        left = root(mass,PLANCK["TARGET"]+PLANCK["HALF_UNIT"],False)
        right = root(mass,PLANCK["TARGET"]-PLANCK["HALF_UNIT"],False)
        dmin = gap(a(left[0])).union(gap(a(left[1])))
        dmax = gap(a(right[0])).union(gap(a(right[1])))
        for sign in (1,-1):
            assert target_gap(gamma,500,sign)<dmin and target_gap(gamma,650,sign)>dmax
            for gain in range(501,650):
                want = target_gap(gamma,gain,sign)
                if want < dmin or want > dmax:
                    continue
                assert dmin<want<dmax
                rlo,rhi = root(gap,want,True)
                r = a(rlo).union(a(rhi))
                scale_fit = mass(r)
                assert abs(scale_fit-PLANCK["TARGET"]) < PLANCK["HALF_UNIT"]
                # Recompute all threshold functions directly at the selected solution.
                s = MATCH["SCATTER"]["scalar_thresholds"](products,r,lam)
                delta = MATCH["threshold_correction"](r*r,sigma)
                low,v,_ = RATIO["MODEL"]["matched"]([v0+dx for v0,dx in zip(s,delta)],eps,k,r)
                assert (v/low*scale-scale_fit).contains(0)
                assert ((a("13/10")*(1+3/k**2)).sqrt()*low*r-want).contains(0)
                yield n,gain,sign,r,scale_fit


def input_sensitivity():
    # The two Yukawa ratios have a common relative shift, so d/h stays fixed.
    e,f = a("1/3072"),a("1/12")
    corner = a(".999").union(a("1.001"))
    p,q = a("7/2")*corner,a("13/10")*corner
    step = a("1/4096")
    margin = 16*(-step).union(step)
    for _ in range(2048):
        if a("16/5")*p.upper()+3*q.upper()<a("206/15") and p.upper()+a("44/5")*q.upper()<a("314/15"):
            break
        u,v = p+margin,q+margin
        du,dv = u*(a("16/5")*u+3*v-a("206/15")),v*(u+a("44/5")*v-a("314/15"))
        assert abs(du)<16 and abs(dv)<16
        p,q = p+step*du,q+step*dv
    else:
        raise AssertionError("No invariant rectangle reached")
    for n,kind,r,v,state in MATCH["matched_points"]():
        if kind!="center":
            continue
        z = 1/state[0]
        want = target_gap(acb.zeta_zero(n).imag,624,1)
        multiplier = (want/(r*(a("13/10")*z).sqrt()))**2
        assert abs(multiplier-1)<a(".000667")
        xi5,xi10 = a("7/2")*multiplier,a("13/10")*multiplier
        a4 = (15-2*xi5**2+6*xi10**2)*z*z+(-60+4*(xi5+3*xi10))*z*e+a("64/5")*e**2-30*f**2
        a22 = (9-6*xi10**2)*z*z+a("168/25")*e**2+30*f**2
        assert a22<0 and a22+a4<0
        print(n,"centered relative ratio shift",(multiplier-1).str(9))


def verify():
    ctx.prec = 150
    results = list(matches())
    input_sensitivity()
    assert [(n,gain,sign) for n,gain,sign,_,_ in results] == [(11490,624,1)]
    for n,gain,sign,r,mass in results:
        orbits = {tuple(sorted(v)) for v in product(range(30),repeat=3) if sum(x*x for x in v)==867}
        words = set(product(orbits,repeat=4))
        words.remove(((17,17,17),)*4)
        assert len(orbits)==5 and len(words)==gain
        print(n,gain,sign,"u/v",r.str(12),"M/GeV",mass.str(12),
              "relative shift",(mass/a(PLANCK["TARGET"])-1).str(8))


if __name__=="__main__":
    verify()
