"""已发表四节点截止在固定带荷谱时唯一识别中性零点高度与有限误差阈值。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx, fmpz_mat


def verify():
    root = Path(__file__).resolve().parent
    common = run_path(str(root/"ch39-five-zero-gravity-gauge-transfer.py"))
    number = common["BASE"]
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    spectral = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))
    single, lib, (rows,) = spectral["gravity_fields"](root, (number,), tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 320
    sites = list(map(arb, common["ORIGINAL_SITES"]))
    gauge = run_path(str(root/"ch39-signed-gravity-unification.py"))
    charges = gauge["gauge_weights"](model)
    columns = [gauge["threshold_column"](rows, single["moments"], lib["expint"], charges, x)
               for x in sites]
    d = (-arb(115)/218, arb(333)/218, arb(-1))
    dot = common["dot"]
    matrix = arb_mat([c[:3]+[dot(d, c[3:])] for c in columns]).transpose()
    weights = matrix.solve(arb_mat([[1], [1], [arb("1e-122")], [-dot(d, initial)]])).entries()
    variation = sum(abs(w) for w in weights)
    assert variation < arb("4.268")

    flat, star = fmpz_mat(19, 19), fmpz_mat(19, 19)
    for i, j, _, _ in shell["geo"]["edges"]:
        tuned = shell["soft"]["pole"] in (i, j) and (i in shell["soft"]["star"]["orbit"]
                                                     or j in shell["soft"]["star"]["orbit"])
        for laplacian in ((flat, star) if tuned else (flat,)):
            laplacian[i, i] += 1
            laplacian[j, j] += 1
            laplacian[i, j] -= 1
            laplacian[j, i] -= 1
    eigen = [(z.real, m) for z, m in flat.charpoly().complex_roots()]
    assert sum(m for v, m in eigen if v.is_zero()) == 1
    identity = arb_mat([[int(i == j) for j in range(19)] for i in range(19)])
    scale = 32*arb("2.435e18")/(724**2*arb("1.22089e19"))
    slope = arb(0)
    for i, (v, mult) in enumerate(eigen):
        if v.is_zero():
            continue
        assert v > 8 and abs(v-15) > arb(".235")
        projector = identity
        for j, (w, _) in enumerate(eigen):
            if j != i:
                projector = projector*(arb_mat(flat)-identity*w)/(v-w)
        product = projector*arb_mat(star)
        drift = -sum(product[j, j] for j in range(19))
        z = (scale*v*(v-15))**2
        slope -= 2*(drift/v+(drift+arb(mult)/3)/(v-15))*sum(
            weight*(-z/x).exp() for weight, x in zip(weights, sites, strict=True))
    assert arb("-1.357") < slope < arb("-1.356")

    # 复参数圆盘中用正规基准矩阵的谱扰动界，保住全部非零模的 Re z>0。
    radius = arb(".01")
    shift = 4*radius
    other = shift+radius/3
    angle = 2*((shift/(8-shift)).atan()+(other/(arb(".235")-other)).atan())
    assert angle < arb(".5")
    lower_real = (scale*(8-shift)*(arb(".235")-other))**2*angle.cos()
    assert lower_real > arb("1e-10")
    bound = 18*variation*arb("1e-10").expint(1)
    offset = shell["soft"]["point"]["study"]["offset"]
    def deficit(height):
        return arb(23)/25*arb(2)**-20/(height**2+offset)
    def derivative_error(loss):
        return bound/radius*((1-loss/radius)**-2-1)
    maximum = deficit(acb.zeta_zero(1).imag)
    assert maximum < arb("3.82e-9") and derivative_error(maximum) < arb(".14")
    assert slope+derivative_error(maximum) < arb("-1.21")
    prev, center, following = [acb.zeta_zero(j).imag for j in (number-1, number, number+1)]
    assert prev < center < following
    central_loss = deficit(center)
    slope_box = slope+arb(0, derivative_error(central_loss).upper())
    projection_limit = -slope_box*central_loss*21/(436*arb.pi())
    epsilon_star = projection_limit/sum(abs(x) for x in d)
    assert arb("1.67539813e-32") < epsilon_star < arb("1.67539815e-32")
    print("中性迹在零亏额的导数、整段导数偏移上界", slope.str(15), derivative_error(maximum).str(10))
    print("基准亏额、高度无界时的投影及逆耦合误差门槛",
          central_loss.str(15), projection_limit.str(15), epsilon_star.str(15))


if __name__ == "__main__":
    verify()
