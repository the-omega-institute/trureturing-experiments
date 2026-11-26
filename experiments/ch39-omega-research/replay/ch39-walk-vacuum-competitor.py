"""原壳两面字公平游走的固定步数及随机停止步数，都有同一参考真空竞争点。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    exact = load("mincut-real-triangle-tomography")
    shell = load("three-ray-chiral-boundary")
    lib = load("mincut-global-winding")
    angles = load("monotone-winding-action")
    ctx.prec = 224
    theta = exact["unit"]((108469,109)).arg()
    bands = [(arb(l),arb(r)) for l,r in angles["BANDS"]]
    dplus, dminus = angles["verify_phases"](root,shell,bands)
    ctx.prec = 224
    tri = lib["terms"](shell)
    cap = arb(23)/25*arb(2)**-20/14**2
    delta = arb(0).union(cap)
    coeffs = [lib["coeff"](tri,q,delta) for q in range(1,106)]
    cosines = [(q*theta).cos() for q in range(1,106)]
    assert cosines[0] > 0
    weights = []
    for n in range(213):
        raw = [z**n for z in cosines]
        normalizer = sum(abs(v) for v in raw)
        assert normalizer > 0
        weights.append([v/normalizer for v in raw])
    assert all(row[-1] < arb(".072") for row in weights[:34:2])
    assert all(row[-1] > arb(".072") for row in weights[34::2])
    print("16、17、106步的末圈份额",*(weights[2*k][-1].str(17) for k in (16,17,106)))

    W = arb_mat(weights)
    # 竞争点两钟都取原生两面字的相角；钟差为零。
    rival = [(acb(0,q*theta).exp(),acb(0,q*theta).exp(),acb(1)) for q in range(1,106)]
    worst_even, worst_odd = arb(100), arb(100)
    leaves = 0
    for sign, d in ((1,dplus),(-1,dminus)):
        for l,r in bands:
            for j in range(64):
                t = (l+(r-l)*j/64).union(l+(r-l)*(j+1)/64)
                for lift in (arb(0),arb.pi()):
                    x, y, diff = sign*t+lift, sign*(t-d)+lift, sign*d
                    values = []
                    for q,(cs,reference) in enumerate(zip(coeffs,rival),1):
                        target = [acb(0,q*s).exp() for s in (x,y,diff)]
                        value = sum((c*(b-a)).real for c,b,a in zip(cs,reference,target))
                        values.append([value])
                    margins = (W*arb_mat(values)).entries()
                    worst_even = min(worst_even, *(v.lower() for v in margins[::2]))
                    worst_odd = min(worst_odd, *(v.lower() for v in margins[1::2]))
                    leaves += 1
    assert worst_even > arb("4.29556") and worst_odd > arb("-2.12130")
    assert (worst_even+worst_odd)/2 > arb("1.08713")
    print("覆盖角盒、偶次试探能差下界、奇次有符号试探能差下界",
          leaves, worst_even.str(17), worst_odd.str(17))
    print("任意停止分布先合成振幅再平方的共同能差下界",
          ((worst_even+worst_odd)/2).str(17))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
