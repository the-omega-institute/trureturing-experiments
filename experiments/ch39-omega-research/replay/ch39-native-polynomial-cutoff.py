"""荷二实际投影和编号二的 Gram 枢轴上，非负系数截止多项式的最小次数。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    power=run_path(str(root/"ch39-native-power-cutoff-matching.py"))
    model=power["native_readout"](root)
    # 在 [0,1]，次数至多 12 的归一非负系数多项式全部满足 r**12 <= P(r) <= 1。
    assert model["floor"]>arb(".009725") and model["low_degree_upper"]<arb("-.06649")
    def nodes(scale,t):
        return [scale*r**12*(1-t*(1-r)) for r in model["ratios"]]

    left,right=arb(".02634403"),arb(".02634404")
    low,high=arb(".99687536"),arb(".99687538")
    scales,ts=left.union(right),low.union(high)
    nleft=model["moments"](nodes(left,ts))[0]-1
    nright=model["moments"](nodes(right,ts))[0]-1
    dlow=model["weak"](nodes(scales,low))
    dhigh=model["weak"](nodes(scales,high))
    print("Newton 和弱差在两对整边上的包围",*[v.str(13) for v in (nleft,nright,dlow,dhigh)],flush=True)
    assert nleft<0<nright and dlow<0<dhigh
    xs=nodes(scales,ts)
    assert all(0<x<1 for x in xs)
    vacuum=model["moments"](xs)[1]
    print("整矩形的 W",vacuum.lower().str(16),vacuum.upper().str(16),flush=True)
    assert arb("-.117246583")<vacuum<arb("-.117246449")
    print("首、末截止质量（GeV）",*[(model["mp"]*x.sqrt()).str(13) for x in (xs[0],xs[-1])])

if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
