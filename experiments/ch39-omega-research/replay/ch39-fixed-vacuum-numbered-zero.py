"""在三枚编号零点上验证固定三角作用的真空能差与味读数误差。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def verify(root, index):
    fixed = run_path(str(root/"ch39-height-independent-flavor-vacuum.py"))
    winding = run_path(str(root/"ch39-mincut-global-winding.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    local = run_path(str(root/"ch39-mincut-reference-ckm-action.py"))
    cp = run_path(str(root/"ch39-three-window-two-probability-cp.py"))["cp"]
    physical = run_path(str(root/"ch39-first-three-mode-observables.py"))
    ctx.prec,ctx.cap = 256,3
    height = acb.zeta_zero(index).imag
    phases = shell["geo"]["phases"](height)
    x,y = [phases[k].arg() for k in (724,804)]
    flavor = cp["flavor"]
    target_x = arb("-.01046")
    target_y = target_x-local["reference_angles"](flavor)[0]
    coordinate = flavor["singlet_at"](phases)[1]
    reference = flavor["singlet"](11661)[1]
    data = cp["two"]["spectral_data"]()
    polynomials = cp["cp_data"](data,flavor["reference_matrix"]())
    def read(c):
        return cp["two"]["probabilities"](c,data)[0][0], cp["cp_value"](c,polynomials)
    readings, benchmark = read(coordinate), read(reference)
    relative = [v/base-1 for v,base in zip(readings,benchmark,strict=True)]
    model,initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    rows = chiral["source"]["spectral_supplier"](model)(height)
    bands = {q:chiral["ordered"]([v for v,k in rows if k == q]) for q in (0,*model["charges"])}
    delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    cut = 9+(1-delta)**2
    assert all(bands[q][2]<9<cut<10<bands[q][3] for q in model["charges"])
    _,alpha,curve = physical["layer_observables"](model,initial,bands)
    assert arb(".117") < alpha < arb(".119") and arb(".998") < curve < arb("1.002")

    tri = winding["terms"](shell)
    def energy_gap(weight_rule):
        weights = weight_rule(winding, tri, target_x, target_y)
        gap = arb(0)
        for q, w in weights.items():
            coefficients = winding["coeff"](tri, q, delta)
            gap += w*(winding["energy"](q, coefficients, x, y)
                      - winding["energy"](q, coefficients, target_x, target_y))
        assert gap > 0
        return gap

    gap = energy_gap(fixed["fixed_weights"])
    optimal_gap = energy_gap(fixed["minimum_order_weights"])
    print("105圈固定作用相对同高度真空能差", optimal_gap.str(17))
    if index == 40_166_104_723:
        assert arb("1.79764315e-5") < optimal_gap < arb("1.79764316e-5")
        assert abs(x-target_x) < arb(".000133") and abs(y-target_y) < arb(".000133")
        assert arb("2.0047874e-8") < coordinate-reference < arb("2.0047875e-8")
        assert arb("-1.313842e-6") < relative[0] < arb("-1.313841e-6")
        assert arb("-8.02236e-7") < relative[1] < arb("-8.02235e-7")
        assert arb("7.723479e-5") < gap < arb("7.723480e-5")
    elif index == 37_179_176_867:
        assert arb("5.79033453e-8") < optimal_gap < arb("5.79033454e-8")
        assert abs(x-target_x) < arb("1e-5") and abs(y-target_y) < arb("1e-5")
        assert abs(coordinate-reference) < arb("9.235e-7")
        assert arb("-6.05203e-5") < relative[0] < arb("-6.05202e-5")
        assert arb("-3.69546e-5") < relative[1] < arb("-3.69545e-5")
        assert arb("1.145990e-7") < gap < arb("1.145991e-7")
    elif index == 2_306_888_100_085:
        # 同一枚编号零点的两个钟相角落入联合密度证书的正支开窗。
        assert abs(y-target_y) < arb("9e-6")
        assert abs(x-y-target_x+target_y) < arb("2e-7")
        assert 0 < optimal_gap < arb("1e-7")
        assert abs(relative[0]) < arb("2e-6")
        assert abs(relative[1]) < arb("1.2e-6")
        print("正支共同开窗内的高度", height.str(30))
    print("零点序号、高度、两相角",index,*(v.str(18) for v in (height,x,y)))
    print("两钟角差、相对角距差、六面坐标差",
          *[v.str(14) for v in (x-target_x,y-target_y,x-y-target_x+target_y,coordinate-reference)])
    print("强耦合、曲率比和第 3、第 4 质量的跨荷分界",
          *[v.str(14) for v in (alpha,curve,max(bands[q][2].upper() for q in model["charges"]),
                                min(bands[q][3].lower() for q in model["charges"]))])
    print("U 与 J 的绝对差、相对差与同高度全局能差",
          *[v.str(14) for v in [*(v-base for v,base in zip(readings,benchmark)),*relative,gap]])


if __name__=="__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    for index in (40_166_104_723, 37_179_176_867, 2_306_888_100_085):
        verify(root,index)
