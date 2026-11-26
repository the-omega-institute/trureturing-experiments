"""Joint threshold, scalar contact and gauge-flow checks for a weaker adjoint self term."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
SCATTER = run_path(str(ROOT/"ch39-quartic-scalar-scattering.py"))
FLOW = run_path(str(ROOT/"ch39-quartic-split-family-flow.py"))
PLANCK, RATIO, a = SCATTER["PLANCK"], SCATTER["RATIO"], SCATTER["a"]


def threshold_correction(q, sigma):
    color = (1+sigma*q)/(1+q)
    return [arb(0), -sigma.log()/6, -color.log()/4]


def charged_eigenvectors():
    import numpy as np
    # Explicit block decoupling on each off-diagonal unbroken adjoint root.
    model = RATIO["MODEL"]
    projector = model["P"]
    color_projector = model["COLOR"]
    end = model["lift"](color_projector)-2*np.eye(10,dtype=int)
    for sign in (1,-1):
        for i in range(5):
            for j in range(i+1,5):
                k = np.zeros((5,5),int)
                k[i,j],k[j,i]=1,sign
                lock = color_projector@k+k@color_projector-k
                if (i<3)!=(j<3):
                    assert not lock.any()
                    continue
                lifted = projector@model["lift"](k)
                assert not (lifted@end.T).any()
                assert (lock**2).sum()==(k**2).sum()
                assert (lifted**2).sum()==((k**2).sum() if j<3 else 0)


def matched_points():
    lam, sigma = a("1/12"),a("1/256")
    products = RATIO["product_polynomials"]()
    powers = PLANCK["mass_exponents"](products)
    match = RATIO["MODEL"]["matched"]
    k = a("8/5")
    qbox = a(".16").union(a(".17"))**2
    slope = ((1344*qbox**3+491*qbox**2-1400*qbox-654)/(42*qbox*(1+qbox)*(2+qbox)*(1+2*qbox))
             -a("107/42")*(sigma/(1+sigma*qbox)-1/(1+qbox)))
    assert slope < 0
    for n,h,eps in RATIO["high_zero_readings"]():
        low,v,_ = match(SCATTER["scalar_thresholds"](products,arb(1),lam),eps,arb(1))
        geom = (2*arb.pi()*(h*h+a("1/4")+30+a("2/15"))/(3*arb(5).sqrt()/2)).sqrt()
        anchor = v/low*geom
        def mass(r):
            q = r*r
            color = (1+sigma*q)/(1+q)
            shift = a("-107/42")*color.log()+a("18/7")*sigma.log()
            return anchor*PLANCK["mass_ratio"](r,powers)*shift.exp()
        def root(target):
            return PLANCK["enclose_root"](mass,target,F(".16"),F(".17"),False)
        left,right = root(PLANCK["TARGET"])
        points = [("center",a(left).union(a(right)))]
        if n==11485:
            fraction = PLANCK["least_denominator"](root(PLANCK["TARGET"]+PLANCK["HALF_UNIT"]),
                                                   root(PLANCK["TARGET"]-PLANCK["HALF_UNIT"]))
            print("Rounded first zero, rational u/v",fraction)
            points.append(("rounded",a(fraction)))
        for kind,r in points:
            base = SCATTER["scalar_thresholds"](products,r,lam)
            delta = threshold_correction(r*r,sigma)
            l,v,lx = match([s+ds for s,ds in zip(base,delta)],eps,k,r)
            direct = v/l*geom
            assert (direct-mass(r)).contains(0)
            if kind=="center":
                assert direct.contains(arb(PLANCK["TARGET"]))
            else:
                assert abs(direct-PLANCK["TARGET"]) < PLANCK["HALF_UNIT"]
            x = lx.exp()
            assert x*eps*v > 1000 and a("1.11") < x < a("1.14")
            state = [1/(l*l*(1+3/k**2)),1/(l*l*(3+k*k)),*([x*x]*3),*([(eps*x)**2]*3)]
            yield n,kind,r,v,state


def verify():
    ctx.prec = 160
    charged_eigenvectors()
    for n,kind,r,v,state in matched_points():
        alpha5,alpha10 = FLOW["full_flow"](state,(a(PLANCK["TARGET"])/v).log())
        assert a(".227") < alpha5 < a(".231") and a(".0179") < alpha10 < a(".0181")
        print(n,kind,"u/v",r.str(9),"v/GeV",v.str(8),"Planck",alpha5.str(7),alpha10.str(7))


if __name__ == "__main__":
    verify()
