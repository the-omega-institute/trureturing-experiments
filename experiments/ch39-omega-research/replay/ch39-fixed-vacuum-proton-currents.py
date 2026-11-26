"""在 105 圈固定真空的联合开窗上验证轻重子电流的 Wilson 因子及其最优下界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


CLASSES = ((1, 0), (0, 1), (1, -1))


def verify(root):
    fixed = run_path(str(root / "ch39-height-independent-flavor-vacuum.py"))
    winding = run_path(str(root / "ch39-mincut-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    angles = run_path(str(root / "ch39-mincut-reference-ckm-action.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    grassmann = run_path(str(root / "ch39-baryon-triality-operator-floor.py"))
    ctx.prec = 224
    diff = angles["reference_angles"](flavor)[0]
    x = arb("-.01046")
    y = x-diff
    tri = winding["terms"](shell)
    groups = [tri[key] for key in CLASSES]
    assert [len(faces) for faces in groups] == [42, 42, 12]
    assert [sum(soft for _, soft in faces) for faces in groups] == [3, 3, 0]
    weights = fixed["minimum_order_weights"](winding, tri, x, y)
    assert max(weights) == 105
    columns = {q: [arb(1), *fixed["gradient_columns"](winding, tri, q, x, y)] for q in weights}
    # 同一解的归一化与四个驻点方程有独立的五行。
    minor = arb_mat([columns[q] for q in (30, 35, 38, 40, 105)])
    assert not minor.det().contains(0)
    for j in range(5):
        total = sum(weight*columns[q][j] for q, weight in weights.items())
        assert (total-int(j == 0)).contains(0) and total.rad() < arb("1e-35")

    def factor(k, angle, delta):
        value, slope = acb(0), acb(0)
        for q, weight in weights.items():
            coefficient = winding["coeff"](tri, q, delta)[k]
            term = weight*coefficient*acb(0, q*angle).exp()/len(groups[k])
            value += term
            slope += acb(0, q)*term
        return value, slope

    cap = arb(23)/25*arb(2)**-20/14**2
    delta = arb(0).union(cap)
    s, t = arb("9e-6"), arb("2e-7")
    window = (x+arb(0,s+t), y+arb(0,s), diff+arb(0,t))
    checked = [factor(k, angle, delta) for k, angle in enumerate(window)]
    bounds = [value.real for value, _ in checked]
    assert bounds[0] > arb(".0801") and bounds[1] > arb(".0555")
    assert arb(".0314") < bounds[2] < arb(".0315")
    assert checked[2][1].real > arb(".969")
    infimum = factor(2, diff-t, arb(0))[0].real
    center = factor(2, diff, arb(0))[0].real
    assert arb(".03143977991031265658") < infimum < arb(".03143977991031265659")
    assert arb(".03143997387646487292") < center < arb(".03143997387646487294")

    # 同一枚编号零点的混合类不含软辐边，实部仅依赖此零点的钟角差。
    height = acb.zeta_zero(2_306_888_100_085).imag
    phases = shell["geo"]["phases"](height)
    actual_diff = phases[724].arg()-phases[804].arg()
    assert abs(actual_diff-diff) < t and abs(phases[804].arg()-y) < s
    numbered = factor(2, actual_diff, arb(0))[0].real
    assert infimum < numbered < center

    # 一代 QQQL 有显式非零 Grassmann 系数。X=2Y-B 模三下 D=QQ 的荷为 0，T=(QL)^* 的荷为 1。
    assert grassmann["qqql"]() == 4
    charge_q, charge_l = 0, -1
    assert (2*charge_q) % 3 == 0 and (-charge_q-charge_l) % 3 == 1
    # 两个等质量辅助场的源向量具有相同范数，混合项只区别于内积。
    for lepton, overlap in [((0, 1), 0), ((1, 0), 1)]:
        diquark = (1, 0)
        assert sum(a*a for a in diquark) == sum(a*a for a in lepton) == 1
        assert sum(a*b for a, b in zip(diquark, lepton)) == overlap
    print("三类三角形在联合开窗上的实部", *(bound.str(11) for bound in bounds))
    print("十二面类的实部导数", checked[2][1].real.str(11))
    print("正反向凸混合的模下确界与真空处最小值", infimum.str(22), center.str(22))
    print("第 2306888100085 枚正零点的最小模", numbered.str(22))
    print("归一化与四个驻点分量满秩，QQQL 单项系数", grassmann["qqql"]())


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
