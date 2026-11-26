"""前105圈共用同一全局真空时，两个指定钟差上的质子因子最小最大模。"""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, arb_mat, ctx


def verify(*, joint_window=False):
    root = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    fixed = load("height-independent-flavor-vacuum")
    model = load("reference-flavor-global-winding")
    lib = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    reference = load("mincut-reference-ckm-action")
    ctx.prec = 224
    tri = lib["terms"](shell)
    model["conjugation_symmetry"](root)
    assert len(tri[1,-1]) == 12 and all(s == 0 for _, s in tri[1,-1])
    d = reference["reference_angles"](load("chiral-neutral-flavor-phase"))[0]
    x = arb("-.01046")
    y = x-d
    number = 2_306_888_100_085
    height = acb.zeta_zero(number).imag
    phases = shell["geo"]["phases"](height)
    dn = (phases[724]/phases[804]).arg()
    radius = arb("2e-7")
    left, right = (d-radius, d+radius) if joint_window else (dn, d)
    gap = d-dn
    assert gap > 0
    # 圆弧各面都是单位相位，因而归一正权下 |F''| <= 105²。
    # 与斜率对偶证书合用：F(d)-F(dn) >= gap*(a-b*F(d))-gap²*105²/2。
    load("proton-clock-difference-slope")["verify"](root)
    a, b = arb(".35511852"), arb("1.7823067441")
    curvature = max(q*q*abs(lib["coeff"](tri, q, arb(0))[2])/12 for q in range(1,106))
    assert curvature < arb("6736.888")
    lower = gap*(a-gap*curvature/2)/(2+b*gap)
    if joint_window:
        # 两端的中心差分消去常数项，截断半宽至多 K*r²/2。
        lower = radius*(a-radius*curvature/2)/(1+b*radius)
        threshold = arb("5e-8")
        slope = a-b*(threshold+105*radius)-curvature*radius
        missed = 1-threshold/(radius*slope)
        assert slope > arb(".3537") and missed > arb(".293")
        print("任意合格权的全窗残差下界、超阈值面积比例下界",
              lower.str(16), missed.str(16))
    else:
        assert lower > arb("2.41e-8")

    # 中点目标行是两个真实读数的算术平均，此处没有用角度中点替代。
    weights = {1: arb(".569438339745"), 32: arb(".003587432201"),
               67: arb(".008363279777"), 69: arb(".001726175641")}
    solved = (34, 35, 38, 83, 104, 105)
    columns, values = {}, {}
    for q in (*weights, *solved):
        coeff = lib["coeff"](tri, q, arb(0))[2]/12
        v, n = [(coeff*acb(0, q*t).exp()).real for t in (right, left)]
        values[q] = (v, n)
        columns[q] = [arb(1), *fixed["gradient_columns"](lib,tri,q,x,y), (v+n)/2]
    mat = arb_mat([columns[q] for q in solved]).transpose()
    assert not mat.det().contains(0)
    rhs = [int(j == 0)-sum(w*columns[q][j] for q, w in weights.items()) for j in range(6)]
    weights.update(zip(solved,mat.solve(arb_mat([[z] for z in rhs])).entries(),strict=True))
    assert min(weights.values()) > arb(".0015")
    delta = arb(0).union(arb(23)/25*arb(2)**-20/14**2)
    modes = [(q,[w*z for z in lib["coeff"](tri,q,delta)]) for q,w in weights.items()]
    drift = model["global_cover"](lib,tri,weights,modes,x,y)
    assert drift < 0
    residual, opposite = [sum(w*values[q][i] for q,w in weights.items()) for i in (0,1)]
    assert (residual+opposite).contains(0)
    assert lower < residual
    if joint_window:
        angle = left.union(right)
        derivative = sum(w*(acb(0, q)*lib["coeff"](tri,q,arb(0))[2]
                            * acb(0,q*angle).exp()).real/12 for q,w in weights.items())
        assert derivative > 0 and residual < arb("7.20e-8")
        # 同一十权作用在整个联合开窗的能差。目标处的梯度精确为零。
        common_shift = arb("9e-6")
        xx = x+arb(0,common_shift+radius)
        yy = y+arb(0,common_shift)
        _, _, hxx, hxy, hyy = lib["derivatives"](modes, xx, yy)
        cap = (abs(hxx+2*hxy+hyy)*common_shift**2 +
               2*abs(hxx+hxy)*common_shift*radius + abs(hxx)*radius**2)/2
        print("十权作用在整个原联合角窗的能差上界",cap.upper().str(15))
        assert cap < arb("1e-7")
    else:
        assert residual < arb("2.46e-8")
    actual_delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    energy = arb(0)
    for q,w in weights.items():
        cs = lib["coeff"](tri,q,actual_delta)
        energy += w*(lib["energy"](q,cs,phases[724].arg(),phases[804].arg())
                     - lib["energy"](q,cs,x,y))
    print("最大残差下界、给定十项正权的残差、实际零点能差",
          *(z.str(18) for z in (lower,residual,energy)))
    assert 0 < energy < arb("1e-7")


if __name__ == "__main__":
    verify()
