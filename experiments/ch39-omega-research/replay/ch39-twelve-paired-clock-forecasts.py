"""对后续十二个不重叠零点对同时封存三种预测，再核对完整指标。"""
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from runpy import run_path

from flint import acb, ctx


TARGETS = tuple(range(23, 46, 2))
BASE = 2_306_888_100_085
SCALE = 10**7


def compute(task):
    root, target, stage = task

    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    shell = lib("three-ray-chiral-boundary")
    phase = shell["geo"]["phases"]
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)[2]
    read_weights = lib("eighth-ninth-wilson-height")["read_weights"]
    arc = lib("eighth-shortest-arc-pivot-support")
    assert arc["SCALE"] == SCALE
    ctx.prec = 320
    new = BASE + target
    target_height = acb.zeta_zero(new).imag
    target_phases = phase(target_height)
    if stage == "actual":
        source = new
        height = target_height
        choices = (target_phases,)
    else:
        source = new - 1
        height = acb.zeta_zero(source).imag
        previous = phase(height)
        choices = (previous, {724: previous[724], 804: target_phases[804]}, target_phases)
    tables = lib("cross-shell-pivot-support")["projected_tables"](
        lib("fourth-third-pivot-global-kinetic-obstruction"), root, source, None
    )
    results = []
    for phases in choices:
        weights = read_weights(shell, gram, height, phases)
        interval = arc["score"](tables, weights)
        assert interval is not None
        results.extend(map(int, interval))
    return (target, *results)


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("stage", choices=("forecast", "actual"))
    options = cli.parse_args()
    root = Path(__file__).resolve().parent
    with ProcessPoolExecutor(max_workers=3) as pool:
        for row in pool.map(compute, ((root, target, options.stage) for target in TARGETS)):
            print("\t".join(map(str, row)), flush=True)
