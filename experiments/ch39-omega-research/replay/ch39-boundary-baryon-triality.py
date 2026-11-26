"""三边界图上代普适的重子数模三、连续规范提升及新增重场阈值。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

boundary = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
model = boundary["model"]


def check_charges():
    # Every renormalizable Yukawa and the Weinberg operator are allowed,
    # while QQQL is charged. Order Q, Uc, Dc, L, Ec, H.
    allowed = []
    for q,u,d,l,e,h in product(range(3),repeat=6):
        terms = (q+u+h, q+d-h, l+e-h, 2*(l+h))
        if all(x%3 == 0 for x in terms) and (3*q+l)%3:
            allowed.append((q,u,d,l,e,h))
    assert len(allowed)==6
    assert all(l==(-h)%3 and h!=0 for q,u,d,l,e,h in allowed)
    # (colour dimension, weak dimension, integer 6Y, X) with X=2Y-B
    sm = [(3,2,1,0),(3,1,-4,-1),(3,1,2,1),(1,2,-3,-1),(1,1,6,2)]
    heavy = [(1,2,3,1),(1,2,-3,2),(1,1,6,-1),(1,1,-6,-2),
             (1,1,0,-3),(1,1,0,0)]
    def anomalies(fields):
        return (sum(w*x for c,w,y,x in fields if c==3), sum(c*x for c,w,y,x in fields if w==2),
                sum(c*w*y*y*x for c,w,y,x in fields),sum(c*w*y*x*x for c,w,y,x in fields),
                sum(c*w*x**3 for c,w,y,x in fields),sum(c*w*x for c,w,y,x in fields))
    assert anomalies(3*sm+heavy)==(0,0,0,0,0,0)
    assert sum(c for c,w,_,_ in 3*sm+heavy if w==2)%2==0
    for a,b in zip(heavy[::2],heavy[1::2]):
        assert a[2]+b[2]==0 and abs(a[3]+b[3])==3
    assert all((q-3*x-b)==0 for (_,_,q,x),b in zip(sm,(1,-1,-1,0,0)))
    print("可用模三荷数、轻谱和抵消场的连续混合/三次/引力反常",len(allowed),anomalies(3*sm),anomalies(heavy))


def verify():
    ctx.prec = 224
    check_charges()
    height = acb.zeta_zero(11431).imag
    b = boundary["boundary_sites"]()
    levels = {q:boundary["high_squares"](boundary["laplacian"](q,height),b) for q in model["charges"]}
    higgs = boundary["high_squares"](boundary["laplacian"](3,height),[boundary["soft"]["pole"]])
    levels[0] = model["shell"]["shell_spectra"]((0,))[0]
    geo = boundary["geo"]
    radius = geo["theta"](804)/geo["theta"](12)
    mz, pi, beta = model["mz"], arb.pi(), model["beta"]
    lo_log,hi_log = [(arb(x)/mz).log() for x in ("1e13","1e17")]
    y,w,c = model["intercept"](levels,radius,boundary["arc"]["first_arc"](),lo_log,
                               fermion_weight=2,higgs_levels=higgs)
    residual = (115*y-333*w+218*c)/218
    assert residual < arb("-4.96")
    # Even a separate light doublet and arbitrarily heavy charged singlet
    # cannot bridge the deficit, as each finite exposure is <= log(mu/M_Z).
    assert residual + 176*hi_log/(218*2*pi) < 0
    for triplets in (0,1,2):
        # A Dirac (1,2)_{1/2}, a Dirac (1,1)_1 and k Weyl (1,3)_0.
        hyper,weak = arb(6)/5,(2+4*arb(triplets))/3
        exposure = -residual*2*pi*218/(333*weak-115*hyper)
        log = (2*pi*(y-w)+(weak-hyper)*exposure)/(beta[0]-beta[1])
        scale = mz*log.exp()
        if triplets==0:
            assert exposure > hi_log
            continue
        if triplets==1:
            assert scale < arb("6.55e12")
            continue
        assert lo_log < log < hi_log and 0 < exposure < log
        common = y-(beta[0]*log+hyper*exposure)/(2*pi)
        assert common > 0
        # With fixed dimensionless M/mu the subtracted response decreases
        # strictly as M increases. Evaluate both sides with endpoint masses.
        def error(ratio):
            kernel = model["sheet"]
            reference_bound = (-2*log).exp()/(8*ratio**2)
            return kernel(ratio**2,arb(1),False)-arb(0).union(reference_bound.upper())-exposure
        low,high = arb("1e-5"),arb(".01")
        assert error(low)>0>error(high)
        for _ in range(38):
            middle = (low+high)/2
            if error(middle)>0:
                low = middle
            else:
                assert error(middle)<0
                high = middle
        mass = low.union(high)*scale
        assert arb("1.2365e13") < scale < arb("1.2367e13")
        print("共同重质量的最少附加弱三重态数、匹配能标、重质量、共同逆耦合",
              triplets,scale.str(12),mass.str(12),common.str(12))
        return scale, mass


if __name__ == "__main__":
    verify()
