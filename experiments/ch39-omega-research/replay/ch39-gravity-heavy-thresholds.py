"""十六重模型中同一非负截止拟合引力时的重层质量源与规范动能阈值。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


def verify(root):
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    single, lib, (rows,) = fields(root, (2_306_888_100_085,))
    model, *_ = lib["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 512
    pi, e, eps = arb.pi(), arb(1).exp(), arb("1e-122")
    k, w = 1/(12*pi), 31/(4*pi)
    heavy = rows[1:]
    assert len(heavy) == 98 and rows[0][0] == -62 and (rows[0][1]-arb(1)/6).contains(0)
    assert rows[0][2].is_zero()
    floor = min(z.lower() for _, _, z in heavy)
    assert floor > arb("1.61822e-11")

    # 将负自旋提前、正自旋延迟，给真实有序前缀一条统一下界。
    events = sorted(((z.upper() if g < 0 else z.lower(), -g) for g, _, z in heavy))
    prefix = 0
    for _, delta in events:
        prefix += delta
        assert prefix >= 0

    budget = eps/(2*w)
    activation = [budget/(e*z) for _, _, z in heavy]
    logarithmic = [4*budget/(e*z)**2 for _, _, z in heavy]
    weighted = sum(abs(a)*v for (_, a, _), v in zip(heavy, activation, strict=True))
    gravity = weighted/(2*pi)
    fraction = max(v.upper() for v in activation)*k/(1-gravity)
    assert max(v.upper() for v in activation) < arb("4.608e-113")
    assert weighted < arb("4.369e-112") and gravity < arb("6.954e-113")
    assert fraction < arb("1.223e-114")

    # 同一固有时处方下的 a4 系数。荷 3 的每个 Dirac 双重态与一个重标量双重态同质量。
    assert tuple(model["charges"]) == (1, 4, 2, 3, 6)
    gauge_weights = []
    scalar_rows = []
    for q, weyl in zip(model["charges"], model["weyl"], strict=True):
        first = len(gauge_weights)
        weights = [2*b + (h if q == 3 else 0) for b, h in zip(weyl, model["higgs"], strict=True)]
        gauge_weights.extend([weights]*16)
        if q == 3:
            scalar_rows.extend(range(first, first+16))
    gauge_weights.extend([(0, -7, -arb(21)/2)]*18)
    gauge = [sum(abs(b[i])*v for b, v in zip(gauge_weights, logarithmic, strict=True))/(4*pi)
             for i in range(3)]
    assert all(v < arb(cap) for v, cap in zip(gauge, ("3.681e-103", "1.689e-102", "8.067e-103"), strict=True))
    portal = sum(activation[j] for j in scalar_rows)/pi
    assert portal < arb("4.022e-113")

    # 四个正截止。两段未知非负高度解具有精确目标的 2x2 系统。
    t = eps/10000
    s = t/2
    u, v = arb("8.1e-12"), arb("1.62e-11")
    assert 0 < s < t < u < v < 1
    # e^(-z/t) <= t²/z²；在极小端点给单截止 K,W 使用相对误差包围。
    tiny = arb("1e-175")
    multiplier = max(6*sum(abs(a) for _, a, _ in heavy), sum(abs(g) for g, _, _ in heavy)/arb(62))
    assert 2*multiplier*(t/floor)**2 < tiny
    def infrared(x):
        ratio = 1 + arb(0, tiny)
        return k*x*ratio, w*x*x*ratio

    ks, ws = infrared(s)
    kt, wt = infrared(t)
    ku, wu = single["moments"](rows, lib["expint"], u)
    kv, wv = single["moments"](rows, lib["expint"], v)
    scale1, scale2 = 1/(k*(t-s)), eps/(w*(v*v-u*u))
    matrix = arb_mat([[(kt-ks)*scale1, (kv-ku)*scale2],
                      [(wt-ws)*scale1/eps, (wv-wu)*scale2/eps]])
    assert not matrix.det().contains(0)
    a, b = matrix.solve(arb_mat([[1-ks], [(eps-ws)/eps]])).entries()
    a, b = a*scale1, b*scale2
    assert a > 1 and b > 0
    assert all(bi[0] >= 0 for bi in gauge_weights)
    gauge_witness = b*(v/u).log()*sum(bi[0]*(-z/u).exp()
                    for bi, (_, _, z) in zip(gauge_weights, heavy, strict=True))/(4*pi)
    assert gauge_witness > arb("1.210e-103")

    print("最轻重层质量平方比、逐层响应、归一逐层响应、绝对加权和、牛顿偏差",
          floor.str(16), max(a.upper() for a in activation).str(16),
          fraction.str(16), weighted.str(16), gravity.str(16))
    print("十六个重标量双重态的单位质量源上界、三级规范动能绝对阈值上界",
          portal.str(16), *[x.str(16) for x in gauge])
    print("精确拟合见证的两段正高度、超荷绝对阈值下界", a.str(16), b.str(16), gauge_witness.str(16))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
