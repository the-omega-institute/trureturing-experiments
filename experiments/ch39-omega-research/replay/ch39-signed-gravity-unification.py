"""同一有符号普朗克有界截止精确匹配牛顿、真空与三规范一圈交汇。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


def verify(root):
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    single, lib, (rows,) = fields(root, (2_306_888_100_085,))
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 320
    pi, eps = arb.pi(), arb("1e-122")
    beta = [b/(2*pi) for b in model["beta"]]
    ratio = (beta[0]-beta[2])/(beta[0]-beta[1])
    d = (1-ratio, ratio, -1)
    def project(v):
        return sum(c*a for c, a in zip(d, v, strict=True))

    weights = []
    assert tuple(model["charges"]) == (1, 4, 2, 3, 6)
    for q, weyl in zip(model["charges"], model["weyl"], strict=True):
        weights.extend([[2*b+(h if q == 3 else 0) for b, h in zip(weyl, model["higgs"], strict=True)]]*16)
    weights.extend([(0, -7, -arb(21)/2)]*18)
    heat = lib["expint"]
    def column(x):
        k, w = single["moments"](rows, heat, x)
        js = [heat(1, z/x) for _, _, z in rows[1:]]
        shifts = [-sum(b[i]*j for b, j in zip(weights, js, strict=True))/(4*pi) for i in range(3)]
        return [arb(1), k, w, *shifts]

    def unification(cols, p):
        intercept = [v+sum(c[3+i]*u for c, u in zip(cols, p, strict=True)) for i, v in enumerate(initial)]
        logscale = (intercept[0]-intercept[1])/(beta[0]-beta[1])
        return project(intercept), model["mz"]*logscale.exp(), intercept[0]-beta[0]*logscale

    # 非负尾函数的规范障碍，包括任何只有一枚截止的测度。
    prefix = 0
    for _, g in sorted((z.upper() if g < 0 else z.lower(), -g) for g, _, z in rows[1:]):
        prefix += g
        assert prefix >= 0
    budget = 4*pi*eps/62
    cap = sum(abs(project(b))*4*budget/(arb(1).exp()*z)**2
              for b, (_, _, z) in zip(weights, rows[1:], strict=True))/(4*pi)
    assert cap < arb("1e-101") and abs(project(initial)) > 4

    # 若两个正截止 a<b 的右端权为负，由 W 目标直接解出该权。
    cut = arb(".001")
    all_k = sum(abs(a) for _, a, _ in rows)/(2*pi)
    pos_k = sum(a for _, a, _ in rows if a > 0)/(2*pi)
    all_spin = sum(abs(g) for g, _, _ in rows)
    small = cut*(pos_k + all_spin*all_k/(2*62))
    assert small < arb(".425")
    kt = sum(a for _, a, _ in rows)/(2*pi)
    wt = -sum(g for g, _, _ in rows)/(8*pi)
    dk = sum(abs(a)*(1-heat(2, z/cut)) for _, a, z in rows)/(2*pi)
    dw = sum(abs(g)*(arb(1)/2-heat(3, z/cut)) for g, _, z in rows)/(4*pi)
    ek = sum(abs(a)*z for _, a, z in rows)/(2*pi*cut)
    ew = sum(abs(g)*z for g, _, z in rows)/(4*pi*cut)
    assert 2*wt > ew
    large = kt+dk+(wt+dw)*(-kt+ek)/(2*wt-ew)
    assert large < arb("-1.285")

    # 三原子。选窄的有理端点计算相反符号，再对整窗复核可逆性及所有物理输出。
    c1, c3 = column(arb("1e-14")), column(arb("3.6e-6"))
    def three(x):
        cols = [c1, column(x), c3]
        m = arb_mat([c[:3] for c in cols]).transpose()
        assert not m.det().contains(0)
        p = m.solve(arb_mat([[1], [1], [eps]])).entries()
        return (*unification(cols, p), p)
    left, right = arb("6.88522069953e-13"), arb("6.88522069955e-13")
    assert three(left)[0] > 0 > three(right)[0]
    _, scale, common, p = three(left.union(right))
    assert p[0] < 0 < p[1]+p[2] and p[2] < 0
    assert arb("5.55606e13") < -p[0]-p[2] < arb("5.55607e13")
    assert arb("3.42999e16") < scale < arb("3.43000e16") and arb("35.9793") < common < arb("35.9794")
    print("排除非负尾的投影与阈值上界、双原子两段上界", project(initial).str(14), cap.str(14),
          small.str(14), large.str(14))
    print("三原子权、交汇能标和统一逆耦合", *[u.str(15) for u in p], scale.str(16), common.str(16))

    # 四个固定有理节点，直接令剩余规范投影为零。
    sites = list(map(arb, ("6.5e-13", ".03", ".5", "1")))
    cols = [column(x) for x in sites]
    m = arb_mat([c[:3]+[project(c[3:])] for c in cols]).transpose()
    assert not m.det().contains(0)
    p = m.solve(arb_mat([[1], [1], [eps], [-project(initial)]])).entries()
    _, scale, common = unification(cols, p)
    assert p[0] > 0 and p[1] > 0 and p[2] < 0 < p[3]
    assert sum(p[1:]) > 0 > sum(p[2:])
    assert arb("1.63377099") < -p[2] < arb("1.63377101")
    assert arb("5.03052e16") < scale < arb("5.03053e16") and arb("37.8345") < common < arb("37.8346")
    print("四原子权、交汇能标和统一逆耦合", *[u.str(17) for u in p], scale.str(17), common.str(17))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
