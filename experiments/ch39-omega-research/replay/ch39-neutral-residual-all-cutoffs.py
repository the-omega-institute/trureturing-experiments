"""归一有符号截止下，全能标窗的独立弱三分量残差及其取消代价。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    axial = run_path(str(root/"ch39-neutral-axial-top.py"))
    boson = run_path(str(root/"ch39-neutral-bosonic.py"))
    ctx.prec = 256
    g1,g2,yt = map(arb,(".467453",".63811",".8616"))
    v,planck,mz = map(arb,("246.22","1.22089e19","91.1876"))
    r = (v/planck)**2/2
    w,z,t = [r*a for a in (g2**2/2,(g2**2+3*g1**2/5)/2,yt**2)]
    p0 = (mz/planck)**2
    pmin = (arb("1e13")/planck)**2
    pi = arb.pi()
    def other_low(p):
        return (-(2*axial["J"](w/p)+boson["j_neutral"](z/p))/(4*pi)
                +3*axial["J"](t/p)/(8*pi)-axial["L"](t/p)/(24*pi))
    low = other_low(p0)+boson["scalar_difference"](w/p0,z/p0)/(8*pi)
    def entropy_bound(mass):
        a = mass/pmin
        return a*(1+(1/(4*a)).log1p())
    # 各个正端点矩随 z/p 增大。对带号组合分别夹住其高能项。
    high = ((2*axial["J"](w/pmin)+boson["j_neutral"](z/pmin))/(4*pi)
            +3*axial["J"](t/pmin)/(8*pi)+axial["L"](t/pmin)/(24*pi)
            +(entropy_bound(w)+entropy_bound(z))/(8*pi))
    coeff = (2*w+z)/(24*pi)+(w+z)/(48*pi)+t/(16*pi)+t/(24*pi)
    inverse_bound = (-low-high)/coeff
    x_min = arb("4e-15")
    tv = inverse_bound*x_min
    error = high+25*coeff/x_min
    assert arb("-.074256577194663365546") < low < arb("-.074256577194663365544")
    assert high < arb("2.398e-21") and coeff < arb("7.468e-36")
    assert inverse_bound > arb("9.9443061e33") and tv > arb("3.9777224e19")
    assert error < arb("4.907e-20") and 25/inverse_bound < arb("2.515e-33")
    print("第八残差参考、全窗高能余量、逆矩系数、逆矩必要下界、原支撑下限必要总变差、T≤25的误差上界",
          low.str(24),high.upper().str(12),coeff.str(18),inverse_bound.lower().str(15),
          tv.lower().str(14),error.upper().str(10))


if __name__ == "__main__":
    verify()
