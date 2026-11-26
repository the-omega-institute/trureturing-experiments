"""弱中心奇偶两个互逆壳刚度的一圈匹配与独立质量标定。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

directory = Path(__file__).parent
shell = run_path(str(directory/"ch39-unit-seventy-five-crossing.py"))
scalar = run_path(str(directory/"ch39-native-shell-higgs-budget.py"))["smooth"]
sheet = shell["smooth"]["sheet"]
q = run_path(str(directory/"ch39-two-loop-shell-separation.py"))["q"]
charges = (1,4,2,3,6)
weyl = [list(map(q,row.split())) for row in
        ("1/15 1 2/3","8/15 0 1/3","2/15 0 1/3","1/5 1/3 0","2/5 0 0")]
higgs = [q("1/10"),q("1/6"),arb(0)]
beta = [q("41/10"),q("-19/6"),arb(-7)]
mz = arb("91.1876")


def setup(zero_index=11431):
    all_levels = shell["shell_spectra"]((0,*charges), zero_index=zero_index)
    levels = {0:all_levels[0]}
    for charge in charges:
        ordered = sorted(all_levels[charge],key=float)
        assert ordered[0] < ordered[1]
        levels[charge] = ordered[1:]
    geo = shell["soft"]["geo"]
    radius = geo["theta"](804)/geo["theta"](12)
    assert all(len(row)==18 for row in levels.values())
    return levels,radius


def intercept(levels,radius,u,Lmin,factor=1,color_inverse=None,fermion_weight=3,higgs_levels=None):
    em,weak = arb("127.930"),arb(".23129")
    if color_inverse is None:
        color_inverse = 1/arb(".1180")
    initial = (3*(1-weak)*em/5,weak*em,color_inverse)
    ref2 = radius**2*(-2*Lmin).exp()
    def integral(charge,kind,spectrum=None):
        if spectrum is None:
            spectrum = levels[charge]
        masses = [v/u if charge in (1,3) else v*u for v in spectrum]
        kernel = {"H":scalar,"V":lambda v,t:sheet(v,t,True),"F":lambda v,t:sheet(v,t,False)}[kind]
        # All angular measures have unit mass and (1-y²)/4 <= 1/4.
        bound = sum(ref2/(8*v) for v in masses)
        return sum(kernel(v,(radius*factor)**2) for v in masses)-arb(0).union(bound.upper())
    f = [integral(charge,"F") for charge in charges]
    v,h = integral(0,"V"),integral(3,"H",higgs_levels)
    return [a+7*n*v/(4*arb.pi())-(b*h+fermion_weight*sum(row[i]*x for row,x in zip(weyl,f)))/(2*arb.pi())
            for i,(a,n,b) in enumerate(zip(initial,(0,2,3),higgs))]


def verify():
    ctx.prec = 256
    levels,r = setup()
    Lmin,Lmax = ((arb(mu)/mz).log() for mu in ("1e13","1e17"))
    def residual(u):
        y,w,c = intercept(levels,r,u,Lmin)
        return (115*y-333*w+218*c)/218
    lo,hi = arb(1),arb(5)
    assert residual(lo)<0<residual(hi)
    for _ in range(32):
        mid = (lo+hi)/2
        value = residual(mid)
        if value<0:
            lo = mid
        else:
            assert value>0
            hi = mid
    u = lo.union(hi)
    y,w,c = intercept(levels,r,u,Lmin)
    L = 2*arb.pi()*(y-w)/(beta[0]-beta[1])
    common = y-beta[0]*L/(2*arb.pi())
    mu = mz*L.exp()
    assert Lmin < L < Lmax and common > 0
    print("Even stiffness, even/odd stiffness ratio, matching scale and common inverse",
          u.str(12),(u*u).str(12),mu.str(12),common.str(12))
    # A quantitative slope supplied by just one odd quark mode is enough.
    x = levels[1][0]/u
    zmin = 3*r*r/16
    low_score = zmin/(x+zmin)
    ref = r*r*(-2*Lmin).exp()
    reference_fraction = ref/(4*x)
    p_derivative = 3*q("180/218")/(2*arb.pi()) * q("11/32")*(low_score-reference_fraction)
    assert p_derivative > arb(".04")
    min_mass = min((v/u if key in (1,3) else v*u).lower()
                   for key,row in levels.items() for v in row)
    ref_fraction = ref/(4*min_mass)
    # 18 modes for each term, coefficients bounded from the tabulated rows.
    coefficient_bounds = [18*(7*n/4+(abs(higgs[i])+3*sum(row[i] for row in weyl))/2)/arb.pi()
                          for i,n in enumerate((0,2,3))]
    assert max(coefficient_bounds) < 45
    assert 180*ref_fraction < arb("1e-17")
    assert (beta[0]-beta[1])/(2*arb.pi())>1
    assert 90*arb("1e-17")/p_derivative < arb("1e-12")
    assert arb("2.83087012") < u < arb("2.83087013")
    assert arb("5.99055e13") < mu < arb("5.99056e13")
    assert arb("38.15697") < common < arb("38.15698")
    def hypercharge(factor):
        y_high = intercept(levels,r,u,Lmin,factor)[0]
        return y_high-beta[0]*(L+factor.log())/(2*arb.pi())
    left,right = arb(1),arb(100)
    assert hypercharge(left)>0>hypercharge(right)
    for _ in range(22):
        mid = (left+right)/2
        value = hypercharge(mid)
        if value>0:
            left = mid
        else:
            assert value<0
            right = mid
    print("Hypercharge inverse first zero in units of matching scale",left.str(11),right.str(11))
    return u, mu, common


if __name__=="__main__":
    verify()
