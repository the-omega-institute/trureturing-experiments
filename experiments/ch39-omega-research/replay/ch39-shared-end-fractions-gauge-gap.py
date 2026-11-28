"""两组实际壳编号及重谱共用逐线端点份额时，对两处规范交汇的分离界。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


def verify():
    def load(name):
        return run_path(str(ROOT / f"ch39-{name}.py"))

    shell = load("native-rh-ray-pvm")["shell_weights"]
    heavy_codes = (2_345_666, 2_345_668)
    _, lib, backgrounds = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT, heavy_codes, tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    heat = lib["expint"]
    weights = load("signed-gravity-unification")["gauge_weights"](model)
    moments = load("signed-cutoff-gravity")["moments"]
    gravity_rows = []
    gravity_center = []
    direction = [arb(n)/100000 for n in (-29071, 84223, 34511, -100000)]
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    slopes = (beta[0]-beta[2], beta[1]-beta[2])
    ref = arb(10**16)
    left, right = (model["mz"]/ref).log(), (arb("1.22089e19")/ref).log()
    center = arb(0)
    worst_shift = arb(0)
    widths = [arb(0) for _ in range(361)]
    for n, (code, rows) in enumerate(zip((2_306_888_100_085, 2_306_888_100_086),
                                         backgrounds, strict=True)):
        y1, y2 = direction[2*n:2*n+2]
        rate = y1*slopes[0]+y2*slopes[1]
        worst_shift += (rate*left).max(rate*right)
        center += y1*(initial[0]-initial[2])+y2*(initial[1]-initial[2])+rate*left
        projected = [y1*(b[0]-b[2])+y2*(b[1]-b[2]) for b in weights]
        pairs = {}
        for charge, levels, direct, _ in shell(ROOT, number=code):
            if charge in (-4, 2):
                pivot = 0 if charge == -4 else 2
                nodes = [v/levels[-1] for v in levels[:-1]]+[arb(1)]
                pairs[charge] = list(zip(nodes, [w.real for w in direct[pivot]], strict=True))
        new_rows = [[], []]
        new_center = [arb(0), arb(0)]
        for i, ((s,w), (t,v)) in enumerate(product(pairs[-4],pairs[2])):
            def response(x):
                return -(w*v)*sum(b*heat(1,z/x) for b,(_,_,z) in
                                   zip(projected,rows[1:],strict=True))/(4*arb.pi())
            x_low = arb(7)/1000*s*t
            x_high = arb(7)/1000*(s+t)/2
            low_pair, high_pair = moments(rows,heat,x_low), moments(rows,heat,x_high)
            for j,(low,high) in enumerate(zip(low_pair,high_pair,strict=True)):
                new_rows[j].append(w*v*(low-high))
                new_center[j] += w*v*(low+high)/2
            first = response(x_low)
            last = response(x_high)
            center += (first+last)/2
            widths[i] += first-last

        gravity_rows.extend(new_rows)
        gravity_center.extend(new_center)

    gravity = arb_mat(gravity_rows)
    gram = gravity*gravity.transpose()
    assert gram.det() > arb("1.937e-12")
    target = arb_mat([[1],[arb("1e-122")],[1],[arb("1e-122")]])
    shift = gravity.transpose()*gram.solve(target-arb_mat([[u] for u in gravity_center]))
    fractions = [arb(1)/2+u for u in shift.entries()]
    assert all(arb(".171") < f < arb(".947") for f in fractions)
    print("共同引力解中所有左端点份额的最小值与最大值",
          min(v.lower() for v in fractions).str(14), max(v.upper() for v in fractions).str(14))
    radius = sum(abs(d) for d in widths)/2
    norm = sum(abs(y) for y in direction)
    fixed = (center-radius)/norm
    free = (center-radius-worst_shift)/norm
    assert fixed > arb(".5398") and free > arb(".5336")
    print("共用 361 个端点份额的规范差最大绝对值下界，固定能标及两组独立选能标",
          fixed.str(14), free.str(14))


if __name__ == "__main__":
    verify()
