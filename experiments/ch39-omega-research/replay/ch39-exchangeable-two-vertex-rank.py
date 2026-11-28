"""六个实际零点排除全部 644,680 对几何顶点之间的连续混合。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
import sys

from flint import arb_mat, ctx


COMPETITORS = ((-4, 111), (2, 19361), (-4, 11035), (-4, 11757), (2, 31330), (-4, 13300))
CENTRAL = {-4: 2345666, 2: 964890}


def integer_gaps(replay):
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    shell = load("three-ray-chiral-boundary")
    pump = load("native-light-projector-pump")
    selector = load("geometric-selector-mixture")
    grouped = load("geometric-exchangeable-priors")
    light = load("light-spectrum-zero-selector")["light_spectra"](replay)
    ctx.prec = 256
    triple = selector["geometry"](shell, pump)
    labels = list(product(grouped["MASKS"], grouped["MODES"]))
    row = {label: i for i, label in enumerate(labels)}
    geometries = grouped["geometry_classes"]()
    factors = []
    for charge in CENTRAL:
        polys = [selector["squared_factors"](triple(charge, *label)) for label in labels]
        factors.append(arb_mat([[p[k] for k in range(7)] for p in polys]))

    def moments(charge, index):
        levels = [energy / 9 - 1 for energy in light(charge, index)]
        return [sum(x**k for x in levels) for k in range(13)]

    central = {charge: moments(charge, n) for charge, n in CENTRAL.items()}
    for charge, n in COMPETITORS:
        delta = [9**12 * (a-b) for a,b in zip(central[charge], moments(charge, n))]
        kernel = arb_mat([[delta[i+j] for j in range(7)] for i in range(7)])
        masks = factors[0] * kernel * factors[1].transpose()
        bounds = []
        for pairs in geometries.values():
            for up, down in grouped["MODE_PAIRS"]:
                gap = sum(masks[row[u, up], row[d, down]] for u,d in pairs) / len(pairs)
                bounds.append(int(gap.lower().floor().unique_fmpz()))
        assert len(bounds) == 1136
        yield bounds


def first_exclusion(left, right, threshold):
    # 尝试同时满足所有仿射差额 ≤ threshold，保留精确有理闭区间。
    low_num, low_den = 0, 1
    high_num, high_den = 1, 1
    for k, (a, b) in enumerate(zip(left, right)):
        x, y = a-threshold, b-threshold
        if min(x, y) > 0:
            return k
        if x > 0 and x*low_den > low_num*(x-y):
            low_num, low_den = x, x-y
        if y > 0 and (-x)*high_den < high_num*(y-x):
            high_num, high_den = -x, y-x
        if low_num*high_den > high_num*low_den:
            return k
    return None


def verify(replay):
    gaps = list(zip(*integer_gaps(replay)))
    removed = [0] * len(COMPETITORS)
    for i, j in combinations(range(len(gaps)), 2):
        first = first_exclusion(gaps[i], gaps[j], 100)
        assert first is not None, (i, j)
        removed[first] += 1
    assert sum(removed) == 644680
    print("按见证次序排除的顶点对数", removed, "统一反转作用差严格大于", 100)


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent)
