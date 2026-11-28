"""在指定上型中央背景，校准六个有符号截止权，并保留 M_Z 处的强逆耦合。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify():
    root = Path(__file__).resolve().parent
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    single, lib, (rows,) = fields(root, (2_345_666,), tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 256
    gauge = run_path(str(root / "ch39-signed-gravity-unification.py"))
    coefficients = gauge["gauge_weights"](model)
    sites = [arb(x) for x in ("1e-11", "1e-9", "1e-7", ".003", ".07", "1")]
    columns = [gauge["threshold_column"](rows, single["moments"], lib["expint"], coefficients, x)
               for x in sites]
    matrix = arb_mat(columns).transpose()
    assert not matrix.det().contains(0)

    beta = [b/(2*arb.pi()) for b in model["beta"]]
    logscale = (arb("1e16")/model["mz"]).log()
    common = initial[2]-beta[2]*logscale
    shifts = [common+b*logscale-a for a, b in zip(initial, beta, strict=True)]
    target = arb_mat([[v] for v in (1, 1, arb("1e-122"), *shifts)])
    weights = matrix.solve(target).entries()
    negative = sum((-w).max(0) for w in weights)
    assert arb("28.0781") < negative < arb("28.0782")
    assert beta[2] < beta[1] < 0 < beta[0] and logscale > 0
    minimum = common+beta[2]*logscale
    assert arb("8.4745") < minimum < arb("8.4746")
    assert arb("44.4911") < common < arb("44.4913")
    print("六截止权", *[w.str(15) for w in weights])
    print("交点逆耦合、全程最小值、负权总量", common.str(15), minimum.str(15), negative.str(15))


if __name__ == "__main__":
    verify()
