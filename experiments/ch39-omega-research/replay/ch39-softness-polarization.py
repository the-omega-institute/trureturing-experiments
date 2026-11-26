"""Rematch the Planck center while varying the adjoint scalar coefficient."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

BOSON = run_path(str(Path(__file__).with_name("ch39-heavy-polarization-background.py")))
MATCH,a = BOSON["MATCH"],BOSON["a"]
PLANCK,RATIO = MATCH["PLANCK"],MATCH["RATIO"]


def center_function(gamma,epsilon,polys,powers,blocks):
    lam,k = a("1/12"),a("8/5")
    geometry = (2*arb.pi()*(gamma**2+a("1/4")+30+a("2/15"))/(3*arb(5).sqrt()/2)).sqrt()
    scalar = lambda r: MATCH["SCATTER"]["scalar_thresholds"](polys,r,lam)
    match = RATIO["MODEL"]["matched"]
    base,v0,_ = match(scalar(arb(1)),epsilon,arb(1))
    anchor = v0/base*geometry
    def radius(relative):
        sigma = (1+a(relative))/256
        def mass(r):
            q = r*r
            shift = -a("107/42")*((1+sigma*q)/(1+q)).log()+a("18/7")*sigma.log()
            return anchor*PLANCK["mass_ratio"](r,powers)*shift.exp()
        lo,hi = PLANCK["enclose_root"](mass,PLANCK["TARGET"],F(".16"),F(".17"),False)
        for _ in range(54):
            mid = (lo+hi)/2
            residual = mass(a(mid))-a(PLANCK["TARGET"])
            assert residual>0 or residual<0
            if residual>0:
                lo = mid
            else:
                hi = mid
        return a(lo).union(a(hi))
    def evaluate(relative,upper=None):
        r = radius(relative)
        sigma = (1+a(relative))/256
        if upper is not None:
            r = r.union(radius(upper))
            sigma = sigma.union((1+a(upper))/256)
        delta = MATCH["threshold_correction"](r*r,sigma)
        low,v,_ = match([s+d for s,d in zip(scalar(r),delta)],epsilon,k,r)
        assert (v/low*geometry-a(PLANCK["TARGET"])).contains(0)
        z,w = low**2*(1+3/k**2),low**2*(3+k**2)
        moments,_ = BOSON["heavy_moments"](r,z,w,gamma,blocks,sigma)
        return moments[1]/moments[0]**2,moments[2]/moments[0]**3
    return evaluate


def verify():
    ctx.prec = 196
    products = RATIO["product_polynomials"]()
    powers = PLANCK["mass_exponents"](products)
    blocks = list(BOSON["charged_blocks"]())
    q = a(".16").union(a(".17"))**2
    window = F(1,40000)
    sigma = (1+a(-window).union(a(window)))/256
    slope = ((1344*q**3+491*q**2-1400*q-654)/(42*q*(1+q)*(2+q)*(1+2*q))
             -a("107/42")*(sigma/(1+sigma*q)-1/(1+q)))
    assert slope<0 and a("18/7")-a("107/42")*sigma*q/(1+sigma*q)>0
    target = a("5.36647005")
    previous = None
    for n,gamma,epsilon in RATIO["high_zero_readings"]():
        evaluate = center_function(gamma,epsilon,products,powers,blocks)
        left,right = evaluate(-window),evaluate(window)
        assert left[0]<target<right[0]
        lo,hi = -window,window
        for _ in range(60):
            mid = (lo+hi)/2
            value,_ = evaluate(mid)
            assert value<target or value>target
            if value<target:
                lo = mid
            else:
                hi = mid
        H,J = evaluate(lo,hi)
        assert H.contains(target)
        if previous is not None:
            assert J-previous>a("4e-8")
        previous = J
        print(n,"relative softness ppm",(a(lo).union(a(hi))*1000000).str(10),
              "third normalized moment",J.str(15))


if __name__ == "__main__":
    verify()
