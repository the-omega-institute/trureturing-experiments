"""在每枚零点的径向刚度上用小整数权重给出八位参考味码的全局真空。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def verify(counts, high_precision=False):
    root = Path(__file__).resolve().parent
    winding = run_path(str(root / "ch39-mincut-global-winding.py"))
    common = run_path(str(root / "ch39-reference-flavor-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec = 224
    common["conjugation_symmetry"](root)
    tri = winding["terms"](shell)
    cap = arb(23)/25*arb(2)**-20/14**2
    assert acb.zeta_zero(1).imag > 14 and shell["soft"]["point"]["study"]["offset"] > 0
    delta = arb(0).union(cap)
    weights = {q: arb(n)/sum(counts.values()) for q, n in counts.items()}
    modes = [(q, [w*z for z in winding["coeff"](tri, q, delta)]) for q, w in weights.items()]
    x, y = arb("-.01"), arb(".139")
    for _ in range(5):
        gx, gy, a, b, d = winding["derivatives"](modes, x, y)
        dx, dy = arb_mat([[a,b],[b,d]]).solve(arb_mat([[gx],[gy]])).entries()
        x, y = (x-dx).mid(), (y-dy).mid()
    drift = common["global_cover"](winding, tri, weights, modes, x, y)
    assert drift < 0
    gx, gy, *_ = winding["derivatives"](modes, x, y)
    # 两个可能真空框的 Hessian 均大于 50 I，且镜像框对正软化量已排除。
    residual_norm = (gx.abs_upper()**2+gy.abs_upper()**2).sqrt().abs_upper()
    # 单位辐边镜像简并时用相同的正角反射包围，六面振幅的奇分量极小。
    error = residual_norm/50
    assert error < arb("1e-7")
    xx, yy = x+arb(0,error), y+arb(0,error)
    def reading(sign):
        return flavor["singlet_at"]({724:acb(0,sign*xx).exp(), 804:acb(0,sign*yy).exp()})[1]

    coordinate = reading(1).union(reading(-1))
    reference = flavor["singlet"](11661)[1]
    assert (256*coordinate).floor().unique_fmpz() == (256*reference).floor().unique_fmpz() == 252
    cp = run_path(str(root / "ch39-cp-ensemble-test.py"))
    ctx.prec, ctx.cap = 384, 3
    data = cp["two"]["spectral_data"]()
    polynomials = cp["cp_data"](data, flavor["reference_matrix"]())
    u_ratio = cp["two"]["probabilities"](coordinate, data)[0][0] / cp["two"]["probabilities"](reference, data)[0][0]
    j_ratio = cp["cp_value"](coordinate, polynomials) / cp["cp_value"](reference, polynomials)
    if high_precision:
        print("高精度整数权的归一读数偏移", (u_ratio-1).str(15), (j_ratio-1).str(15))
        assert abs(u_ratio-1) < arb("4.7e-7") and abs(j_ratio-1) < arb("2.9e-7")
    else:
        assert arb("1.0041") < u_ratio < arb("1.0043") and arb("1.002") < j_ratio < arb("1.003")
    print("实际相位 U、J 与参考的比值", u_ratio.str(14), j_ratio.str(14))
    for sign in (-1, 1):
        common["physical_all_heights"](root, sign*xx, sign*yy, cap)
    print("相角、梯度余量半径、六面坐标", x.str(13), y.str(13), error.str(12), coordinate.str(14))


if __name__ == "__main__":
    verify({1: 1, 2: 1, 30: 2, 32: 1, 33: 15, 36: 6})
    verify({1: 258, 2: 161, 30: 566, 32: 377, 33: 4044, 36: 1493}, high_precision=True)
