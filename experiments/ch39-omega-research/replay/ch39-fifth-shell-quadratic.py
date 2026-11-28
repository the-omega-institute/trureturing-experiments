"""以二阶质量展开及其逐场余项，盲检第五壳的六个重谱。"""
from functools import partial
import json
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
HOLDOUT = tuple(range(2_345_688, 2_345_699, 2))
ADDITIONAL = (-7175983059, -100000000000, 1468137571, -4256584471)


def lower_bound(context, edges, rows, components, order):
    fields, beta, limit, initial, kernel, earlier = context[1:]
    k, h, u, v = [arb(t)/10**11 for t in components]
    direction = [u, v, -u-v]
    norm = abs(u)+abs(v)
    scale = arb(7)/1000
    assert max(z.upper() for _, _, z in rows) < min((scale*a).lower() for a, _, _ in edges)
    coeff = arb_mat([[k, h, *direction]]) * kernel["series"](rows, fields, order=order)
    basis = partial(kernel["basis"], order=order)
    slope = u*(beta[0]-beta[2])+v*(beta[1]-beta[2])
    lower = -k-h*arb("1e-122")+u*(initial[0]-initial[2])+v*(initial[1]-initial[2])
    lower -= (slope*limit).max(0)
    for a, b, w in edges:
        lower += earlier["projected_minimum"](basis, kernel["S"], [
            (scale*a, scale*(b-a), w*coeff)], arb("1e-9"))
    # xE2、x²E3 及 E1 的交错尾，各由本场第一遗漏项控制。
    n = order+1
    weighted = arb(0)
    for (g, alpha, z), charges in zip(rows, fields, strict=True):
        gauge = sum(d*q for d, q in zip(direction, charges, strict=True))
        weighted += z**n*(abs(2*k*alpha)*z/(n+1)+abs(h*g)*z**2/((n+1)*(n+2))+abs(gauge))
    geometry = sum(abs(w)/(scale*a)**n for a, _, w in edges)
    error = weighted*geometry/(4*arb.pi()*n*factorial(n))
    return (lower-error).lower()/norm, error.upper()/norm


def verify():
    load = run_path(str(ROOT/"ch39-four-background-crossing.py"))["load"]
    exchange = load("fourth-fifth-shell-exchange")
    context = load("shared-large-scale-gap")["load"]()
    shell = exchange["FIFTH_SHELL"]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT, load("native-rh-ray-pvm")["shell_weights"], shell)
    _, _, spectra = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT, HOLDOUT, tight_neutral=True)
    fixed = exchange["SEPARATOR"]
    excluded = {2_345_688: ".1857", 2_345_692: ".4312", 2_345_694: ".7889", 2_345_696: ".1808"}
    for number, rows in zip(HOLDOUT, spectra, strict=True):
        lower, error = lower_bound(context, edges, rows, fixed, 2)
        full, _ = lower_bound(context, edges, rows, fixed, 48)
        print(number, "冻结方向二阶下界、余项上界、48阶下界",
              lower.str(12), error.str(8), full.str(12), flush=True)
        assert error < arb(".00413")
        if number in excluded:
            assert min(lower, full) > arb(excluded[number])
        else:
            assert max(lower, full) < 0
        if number == 2_345_690:
            revised, _ = lower_bound(context, edges, rows, ADDITIONAL, 2)
            assert revised > arb(".00424")
            print("该格另外选取的分离方向，二阶下界", revised.str(12), flush=True)
    layout = json.loads((ROOT/"ch39-fifth-shell-quadratic-layout.json").read_text())
    assert layout["heavy"] == HOLDOUT[-1]
    labels = [arb(k)/11000 for row in layout["base"] for k in row]
    extra = [(None, i, arb(k)/11000) for i, k in layout["extra"]]
    load("one-background-specific-line")["solve_layout"](
        [labels], extra, 0, ["11.74"], [shell], [layout["heavy"]], False, reference_gev=10**12)


if __name__ == "__main__":
    verify()
