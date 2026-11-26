"""Uniform finite enclosures for the quartic model's two-loop gauge running."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
TEN_ONE, TEN_SELF, TEN_CROSS = run_path(str(ROOT/"ch39-two-loop-product-beta.py"))["target_inverse_coefficients"]()
PLANCK = run_path(str(ROOT/"ch39-quartic-planck-branches.py"))
RATIO = PLANCK["RATIO"]
a = RATIO["a"]


def initial_data():
    products = RATIO["product_polynomials"]()
    powers = PLANCK["mass_exponents"](products)
    match = RATIO["MODEL"]["matched"]
    readings = list(RATIO["high_zero_readings"]())
    eps = readings[0][2].union(readings[-1][2])
    assert all(eps.contains(e) for _,_,e in readings)
    # Prove a common rational window contains every small-branch rounding interval.
    edges = [a(".452829"), a(".452847")]
    for n,height,e in readings:
        l,v,_ = match(RATIO["thresholds"](products,arb(1)), e, arb(1))
        geom = (2*arb.pi()*(height**2+a("1/4")+30+a("2/15"))/(3*arb(5).sqrt()/2)).sqrt()
        bounds = [v/l*geom*PLANCK["mass_ratio"](r,powers) for r in edges]
        assert bounds[0] > PLANCK["TARGET"]+PLANCK["HALF_UNIT"]
        assert bounds[1] < PLANCK["TARGET"]-PLANCK["HALF_UNIT"]
    for j in range(32):
        r = (edges[0]+(edges[1]-edges[0])*j/32).union(edges[0]+(edges[1]-edges[0])*(j+1)/32)
        scores = RATIO["thresholds"](products,r)
        # X increases with k on the k interval used below.
        l,v,lx = match(scores,eps,a("33/20"),r)
        assert a(".7525") < l < a(".7532") and a("6e14") < v < a("6.1e14")
        assert lx.exp() < a(".027") and (1+eps**2)*(2*lx).exp() < a(".001")
        assert lx.exp()*eps*v > 1000
    return a(".7525").union(a(".7532"))


def gauge_rhs(g,h,y):
    # g and h are the inverse squared gauge couplings; y is |x|^2+|y|^2.
    unit = 1/(8*arb.pi()**2)
    return (unit*(a("10/3")-(a("6796/15")/g+594/h-45*y)*unit/2),
            unit*(TEN_ONE+(TEN_SELF/h-TEN_CROSS/g+15*y)*unit/2))


def exclude_quarter(low):
    pi = arb.pi()
    step = a("1/256")
    unit = 1/(8*pi*pi)
    slope = unit*(TEN_ONE-TEN_CROSS*pi*unit/2)
    assert slope > 0
    kmin = (3/(pi/low**2-1)).sqrt()
    kmax = (pi/low**2-3).sqrt()
    assert kmin > a(".80") and kmax < a("1.60")
    # Under the global alpha <= 1/4 hypothesis, H grows at least linearly.
    # Integrate its upper derivative using this linear lower bound.
    hmin = 1/pi
    for kleft,kright in ((a(".80"),a("1.3")),(a("1.3"),a("1.60"))):
        g = (1/(low**2*(1+3/kright**2))).upper()
        h0 = (1/(low**2*(3+kleft**2))).upper()
        assert g < a(1)
        for j in range(2500):
            t = (j+1)*step
            h = (h0+unit*(TEN_ONE+15*a(".001")*unit/2)*t
                 +unit**2*TEN_SELF/(2*slope)*(1+slope*t/hmin).log())
            dg,_ = gauge_rhs(g,h,a(".001"))
            assert dg < 0
            g = (g+step*dg).upper()
            if g < 1/pi:
                break
        assert g < 1/pi and a("6.1e14")*t.exp() < PLANCK["TARGET"]
        print("Excluded k range",kleft.str(4),kright.str(4),"by t",t.str(7))


def upper_witness(low):
    pi = arb.pi()
    k = a("33/20")
    g = 1/(low**2*(1+3/k**2))
    h = 1/(low**2*(3+k*k))
    step = a("1/512")
    radius = (-step).union(step)
    y = arb(0).union(a(".001"))
    stop = (a(PLANCK["TARGET"])/a("6e14")).log()
    nsteps = int((stop/step).ceil().unique_fmpz())
    maximum = a(".31")
    for _ in range(nsteps):
        gt, ht = g+radius, h+radius
        dg,dh = gauge_rhs(gt,ht,y)
        assert gt > 0 and ht > 0 and dg < 0 < dh and -dg < 1 and dh < 1
        # Both upper bounds on the Yukawa components persist while G < 1.
        assert gt < 1 and a("54/5")/gt > 13*a(".001")
        assert 1/(4*pi*gt) < maximum and 1/(4*pi*ht) < maximum
        # G' grows in both gauge variables. H' grows in G and falls in H.
        g_low,g_high = g.lower(),g.upper()
        h_low,h_high = h.lower(),h.upper()
        dg_lo = gauge_rhs(g_low-step,h_low,arb(0))[0]
        dg_hi = gauge_rhs(g_high,h_high+step,a(".001"))[0]
        dh_lo = gauge_rhs(g_low-step,h_low+step,arb(0))[1]
        dh_hi = gauge_rhs(g_high,h_high,a(".001"))[1]
        # For H the bounding solutions use their own H endpoint (comparison).
        g = (g_low+step*dg_lo).union(g_high+step*dg_hi)
        h = (h_low+step*dh_lo).union(h_high+step*dh_hi)
    print("Full-trajectory upper cap",maximum.str(5),"final G,H",g.str(5),h.str(5))


def verify():
    from fractions import Fraction as F
    ctx.prec = 140
    coeff = run_path(str(ROOT/"ch39-two-loop-product-beta.py"))["coefficients"]
    b, diagonal = coeff(5,[(F(3,2),F(12,5)),(F(9,2),F(18,5))],
                        [(30,F(18,5)),(F(1,2),F(12,5)),(F(5,2),5)])
    assert (b,diagonal)==(F(-10,3),F(6796,15))
    low = initial_data()
    exclude_quarter(low)
    upper_witness(low)


if __name__ == "__main__":
    verify()
