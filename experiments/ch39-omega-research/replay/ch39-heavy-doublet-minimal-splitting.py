"""固定十六层正零点重谱与五个截止时，抵消独立中性残差所需的最小相对分裂。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


ROOT = Path(__file__).resolve().parent


def verify():
    ctx.prec = 224
    rotation = run_path(str(ROOT/"ch39-rotation-portal-top-energy.py"))
    gauge = run_path(str(ROOT/"ch39-all-field-signed-crossing.py"))
    axial = run_path(str(ROOT/"ch39-neutral-axial-top.py"))
    boson = run_path(str(ROOT/"ch39-neutral-bosonic.py"))
    signs, unshifted = rotation["rotation_signs"](ROOT, 2_306_888_100_085)
    mp, mz, v, energy, g1, g2, portal = map(arb, (
        "1.22089e19", "91.1876", "246.22", "1.176397e14", ".467453", ".63811", ".42995401872"))
    pi = arb.pi()
    rho = (v/mp)**2/2
    p0, p1 = (mz/mp)**2, (energy/mp)**2
    masses = [z+sign*portal*rho for z,sign in zip(unshifted,signs,strict=True)]
    assert len(masses) == 16 and min(masses) > 0
    xmin, xmax = arb("4e-15"), arb(1)
    rule = gauge["angular_rule"](24,57)
    log_ratio = (p1/p0).log()

    def heavy_gain(epsilon):
        total = arb(0)
        for b in masses:
            for x,weight in ((xmax,13),(xmin,-12)):
                light, heavy = [gauge["angular_exposure"](1,z,x,p0,p1,log_ratio,rule)
                                for z in (b*(1-epsilon),b*(1+epsilon))]
                total += weight*(light-heavy)
        return total/(48*pi)

    # F 的导数符号在整个允许顶质量区间只切换一次。
    def slope(t):
        a = (1+4*t).sqrt()
        return (33*a-15/a)*(1/a).atanh()-33

    left, right = map(arb, (".0859553180260300229",".0859553180260300231"))
    assert slope(left) > 0 > slope(right)
    optimum = left.union(right)
    assert 0 < optimum*p0/rho < 4*pi
    j, l = axial["J"], axial["L"]
    fmax = (9*j(optimum)-l(optimum))/(24*pi)
    w,z = rho*g2**2/2, rho*(g2**2+3*g1**2/5)/2
    low_boson = -(2*j(w/p0)+boson["j_neutral"](z/p0))/(4*pi)
    low_boson += boson["scalar_difference"](w/p0,z/p0)/(8*pi)

    # 有符号权的有限截止余量及全部允许顶质量的高端量。
    def entropy(mass):
        r = mass/p1
        return r*(1+(1/(4*r)).log1p())
    tmax = 4*pi*rho
    high_boson = ((2*j(w/p1)+boson["j_neutral"](z/p1))/(4*pi)
                  +(entropy(w)+entropy(z))/(8*pi))
    tr = tmax/p1
    high_top = (3*tr*(1/(4*tr)).log1p()/8+tr/4)/pi
    finite = 25*((2*w+z)/(24*pi)+(w+z)/(48*pi)+5*tmax/(48*pi))/xmin
    error = high_boson+high_top+finite
    assert error < arb("6e-19")
    residual = low_boson+fmax+arb(0,error.upper())
    assert residual < 0
    lower,upper = map(arb, (".2577142859",".2577142869"))
    below,above = residual+heavy_gain(lower),residual+heavy_gain(upper)
    assert below < arb("-1.137e-10") and above > arb("1.098e-10")
    print("最优低端顶质量平方比、顶耦合、轻场余量",optimum.str(18),
          (optimum*p0/rho).sqrt().str(12),error.upper().str(8))
    print("两个分裂端点的残差、对应的重轻质量比",below.str(8),above.str(8),
          ((1+lower)/(1-lower)).sqrt().str(12))


if __name__ == "__main__":
    verify()
