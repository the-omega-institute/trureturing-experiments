"""固定四项时间和一个夸克概率，检验任意有限精度构造的严格内点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


root = Path(__file__).parent
two = run_path(str(root / "ch39-two-ckm-probabilities.py"))
model = two["joint"]["model"]


def verify():
    # 同一相对相位的三窗开邻域在极限壳上有严格裕量，包括辐边软化误差。
    model["physical"]["joint"]["verify"](keep_flavor=True)
    ctx.prec, ctx.cap = 768, 3
    shared = model["shared"]
    (e1, e2, e3), reference, lambdas, slopes = shared["flavor_data"]()
    times = shared["clock"]["clock"]["times"]()
    inv = two["invariant"]
    flavor = inv["flavor"]
    axes = [flavor["singlet_at"]({724:z, 804:acb(1)})[1] for z in (acb(1), acb(0,1))]
    assert sum(k*k for k in axes) > reference**2
    data = two["spectral_data"]()
    up, _ = two["probabilities"](reference, data)
    assert up[1] < 0

    scale = arb.union(*shared["determinant"]["scale_interval"]())
    nodes = (1, 2529893, 2658737, 4035321, 9896608, 25215595)
    ctx.prec = 224
    heights = [acb.zeta_zero(n).imag for n in nodes]
    ctx.prec = 768
    xs = [(scale/h)**2 for h in heights]
    assert e1.degree() == 1 and e2[2] < 0
    assert e2.derivative().derivative()(arb(0, 1))/xs[0] + e3.derivative().derivative()(arb(0, 1)) < 0
    columns = [[arb(1), *[1/(1+j*j*x) for j in range(1,5)], model["gradient"](x,lambdas,slopes)] for x in xs]
    base = arb_mat(list(zip(*columns, strict=True)))
    weights = base.solve(arb_mat([[v] for v in [*times[:5],0]])).entries()
    assert all(w > 0 for w in weights) and weights[4] > arb(".64")
    # 在理想节点上把第五列拆成相反微扰，小幅度下两个系数均为正。
    shift = arb(2)**-22
    left, _ = two["probabilities"](reference-shift, data)
    right, _ = two["probabilities"](reference+shift, data)
    split = weights[4]*(up[0]-right[0])/(left[0]-right[0])
    assert 0 < split < weights[4]
    print("参考点 U 导数、待拆节点权重、一次拆分的两个权重",
          up[1].str(12),weights[4].str(12),split.str(12),(weights[4]-split).str(12))


if __name__ == "__main__":
    verify()
