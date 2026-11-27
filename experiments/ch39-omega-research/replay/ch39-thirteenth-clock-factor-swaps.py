"""第十二到十三壳，逐项交换两种钟边、响应表和软辐边电导。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    numbers = (2_306_888_100_097, 2_306_888_100_098)
    shell = lib("three-ray-chiral-boundary")
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    support = lib("cross-shell-pivot-support")
    arc = lib("eighth-shortest-arc-pivot-support")
    read = lib("eighth-ninth-wilson-height")["read_weights"]
    tables = [support["projected_tables"](core, root, n, None) for n in numbers]
    ctx.prec = 320
    heights = [acb.zeta_zero(n).imag for n in numbers]
    clocks = [shell["geo"]["phases"](h) for h in heights]
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)[2]
    data = {}
    for a, b, soft in product(range(2), repeat=3):
        phases = {724: clocks[a][724], 804: clocks[b][804]}
        weights = read(shell, gram, heights[soft], phases)
        for response, table in enumerate(tables):
            bounds = arc["score"](table, weights)
            assert bounds is not None
            lo, hi = map(int, bounds)
            assert lo > 0 or hi < 0
            data[a, b, response, soft] = (lo, hi)
            print("724、804、响应表、软辐边，及两界分子", a, b, response, soft, lo, hi, flush=True)
    assert arc["SCALE"] == 10**7 and len(data) == 16
    assert data[0, 1, 0, 0] == (-7967044, -7967040)
    assert data[1, 1, 1, 1] == (13283, 13287)
    for (a, b, response, soft), (lo, hi) in data.items():
        assert (hi < 0) == (a == 0 and b == 1)
        assert (lo, hi) == data[a, b, response, 1-soft]
    for response, soft in product(range(2), repeat=2):
        old, both = data[0, 0, response, soft], data[1, 1, response, soft]
        only724, only804 = data[1, 0, response, soft], data[0, 1, response, soft]
        lower = both[0] - only724[1] - only804[1] + old[0]
        upper = both[1] - only724[0] - only804[0] + old[1]
        assert 7648000 < lower <= upper < 7690000
    return data


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
