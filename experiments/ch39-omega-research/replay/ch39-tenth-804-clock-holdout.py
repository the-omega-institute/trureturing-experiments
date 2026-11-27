"""在第十壳勾股点对照只更新 804 时钟的预测与完整读数。"""
from argparse import ArgumentParser
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def evaluate(root, stage):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    old, new = 2_306_888_100_094, 2_306_888_100_095
    shell = lib("three-ray-chiral-boundary")
    clocks = shell["geo"]["phases"]
    ctx.prec = 320
    h0, h1 = [acb.zeta_zero(n).imag for n in (old, new)]
    if stage == "forecast":
        number, height = old, h0
        phases = {724: clocks(h0)[724], 804: clocks(h1)[804]}
    else:
        number, height, phases = new, h1, clocks(h1)
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)[2]
    weights = lib("eighth-ninth-wilson-height")["read_weights"](shell, gram, height, phases)
    support = lib("cross-shell-pivot-support")
    arc = lib("eighth-shortest-arc-pivot-support")
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    tables = support["projected_tables"](core, root, number, None)
    score = arc["score"](tables, weights)
    assert score is not None
    lo, hi = map(int, score)
    print(stage, "S3 的区间分子、分母", lo, hi, arc["SCALE"], flush=True)
    return arb(lo).union(arb(hi)) / arc["SCALE"]


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("stage", choices=("forecast", "actual"))
    evaluate(Path(__file__).resolve().parent, cli.parse_args().stage)
