"""105圈固定全局真空的十二面质子因子及其钟差斜率的有理对偶证书。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


# 每个整数代表除以 10^10 的有理数。前五项乘归一与四条梯度等式，
# 第六项乘 F。每个三元组为竞争点的两角及其非负能差乘子。
scale = 10**10
certificates = (
    (1, (3551185381,-1173097465,-7030824,-428883011,64558196,-17823067441), (
        (-144186871,353916708,5517499098), (-153923116,-971885260,5869836338),
        (-143437520,352939831,5350396218))),
    (-1, (-12190764358,-759107198,-548166491,-270328968,107169570,-212935810542), (
        (283882871,-26326973674,200285034), (-117087676,-1483795602,23812137),
        (-257557605,492131379,3551503045), (-1638676905,-1360062004,8561008877),
        (-228403090,-17261980241,3646778), (-116923601,-1484161712,65757155939),
        (-257734469,492552280,6893981575), (-210816807,-8189155978,29343792))),
)


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    lib = load("mincut-global-winding")
    fixed = load("height-independent-flavor-vacuum")
    reference = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    ctx.prec = 224
    tri = lib["terms"](load("three-ray-chiral-boundary"))
    assert len(tri[1,-1]) == 12 and all(j == 0 for _,j in tri[1,-1])
    x = arb("-.01046")
    d = reference["reference_angles"](flavor)[0]
    y = x-d
    assert d<0
    for sign, eq, raw_cuts in certificates:
        h = [arb(i)/scale for i in eq]
        cuts = [[arb(v)/scale for v in row] for row in raw_cuts]
        assert all(lam>0 for _,_,lam in cuts)
        bounds = []
        for q in range(1,106):
            coeff = lib["coeff"](tri,q,arb(0))
            z = coeff[2]*acb(0,q*d).exp()/12
            factor, slope = z.real,-q*z.imag
            col = [*fixed["gradient_columns"](lib,tri,q,x,y),factor]
            base = lib["energy"](q,coeff,x,y)
            constraint = sum(lam*(lib["energy"](q,coeff,a,b)-base) for a,b,lam in cuts)
            remaining = sign*slope-sum(a*b for a,b in zip(h[1:],col,strict=True))-constraint
            bounds.append(remaining)
        lower = min(v.lower() for v in bounds)
        assert lower > h[0]-arb(".000001")
        print("斜率符号、竞争点数、仿射截距下界、F 的系数",sign,len(cuts),
              lower.str(16),h[-1].str(16))
        assert lower > arb(".35511852" if sign==1 else "-1.21907654")

    # 当参考因子为零，同一组权在实际见证角的读数。十二面的模均为一。
    height = acb.zeta_zero(2306888100085).imag
    phases = load("arc-relative-determinant")["phases"](height)
    step = (phases[724]/phases[804]).arg()-d
    curvature = max(q*q*abs(lib["coeff"](tri,q,arb(0))[2])/12 for q in range(1,106))
    assert step < 0 and curvature < arb("6736.888")
    assert arb(".35511852")-curvature*abs(step)>arb(".35420")
    relative = arb(".35511852").union(arb("1.21907654"))*step + arb(0,curvature*step**2/2)
    print("实际见证的角差偏移、曲率绝对值界和零因子的实际偏移",step.str(13),curvature.str(13),
          relative.lower().str(13),relative.upper().str(13))
    assert relative > arb("-1.6598650e-7") and relative < arb("-4.8271494e-8")


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
