"""四轮十二个门以有限高斯有理门逼近时，三壳九个读数的共同误差证书。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def unit_pair(numerator, bits):
    denominator = 1 << bits
    total = denominator**2 + numerator**2
    return Fraction(denominator**2 - numerator**2, total), Fraction(2*denominator*numerator, total)


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    gate_data = lib("twelve-independent-gates-nine")
    pairs = gate_data["PAIRS"]
    assert len(pairs)==12 and all(len({i for e in pairs[k:k+3] for i in e})==6 for k in range(0,12,3))
    _,squared = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=2_306_888_100_086)[2]
    low,top = squared[0].sqrt(),squared[-1].sqrt()
    native = lib("native-rh-ray-pvm")
    bits = 449
    for m in (0,1,-1,(1<<(bits-1))-1):
        x,y = unit_pair(m,bits)
        assert x*x+y*y==1

    for number in range(2_306_888_100_092,2_306_888_100_095):
        _,_,pivots,_ = next(row for row in native["shell_weights"](root,number) if row[0]==6)
        variation = sum(abs(w.real) for w in pivots[2])
        _,_,(field,) = lib("three-index-free-atoms")["broken_fields"](root,[number])
        ctx.prec = 256
        assert all(mass >= 0 for _,_,mass in field)
        # 单能端的带号弱积分。对 log x 求导时，每项 E1 变为 0 至 1 的指数。
        rho = (arb("246.22")/arb("91.1876"))**2/2
        w = arb(".63811")**2*rho/2
        z = w+3*arb(".467453")**2*rho/10
        t = arb(".8616")**2*rho
        one_endpoint = (w/2+z/4+arb(1)/32+3*t/8+arb(1)/16)/arb.pi()
        # 两个能端作差，且 x=a²。沿酉路径 |(log a)'|≤sqrt((top-low)/low) ||U'||。
        weak = 16*variation*4*one_endpoint*((top-low)/low).sqrt()
        newton = 32*(top-low)*variation*top*sum(abs(alpha) for _,alpha,_ in field)/arb.pi()
        vacuum = 32*(top-low)*variation*top**3*sum(abs(g) for g,_,_ in field)/(2*arb.pi())
        bound = max(weak.upper(),newton.upper(),vacuum.upper())
        print("壳与半角网格步长前的牛顿、真空、弱差系数",number,
              newton.upper().str(12),vacuum.upper().str(12),weak.upper().str(15),flush=True)
        assert bound < arb("1.344e13") and 1344*10**10*10**122 < 2**bits


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
