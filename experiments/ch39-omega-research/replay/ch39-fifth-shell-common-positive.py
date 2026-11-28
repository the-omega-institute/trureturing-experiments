"""固定第五壳，检验三个分别有正解的重谱能否沿用共同谱线分布。"""
import json
from itertools import combinations
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
HEAVY = (2_345_670, 2_345_676, 2_345_684)
# (D_a, K_a-1, D_a-D_b, K_a-K_b, W_a-W_b) 的整数或有限小数系数。
DIRECTIONS = (
    ("-.0031088", "-.02749002", "1", "51975.50479306", "-10000000"),
    ("-.00047913", "-.09635324", "-1", "1961.2735288", "-10000000"),
    ("-.00314666", "-.06886322", "-1", "19007.65331563", "-10000000"),
)


def verify():
    four = run_path(str(ROOT / "ch39-four-background-crossing.py"))
    load = four["load"]
    ctx = load("shared-large-scale-gap")["load"]()
    fields, beta, initial, kernel = ctx[1], ctx[2], ctx[4], ctx[-2]
    shell = load("fourth-fifth-shell-exchange")["FIFTH_SHELL"]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT, load("native-rh-ray-pvm")["shell_weights"], shell)
    scale = arb(7)/1000
    moments = load("fifth-shell-transfer")["shell_moments"](edges, kernel, scale)
    floor = scale * min(a.lower() for a, _, _ in edges)
    _, _, rows = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT, HEAVY, tight_neutral=True)
    assert all(max(z.upper() for _, _, z in bg) < min(arb(1), floor) for bg in rows)
    series = [kernel["series"](bg, fields) for bg in rows]
    t = (beta[0]-beta[2])/(beta[0]-beta[1])
    intercept = (1-t)*initial[0] + t*initial[1] - initial[2]
    gauge_norm = abs(1-t) + abs(t)
    prefactor = 2*sum(abs(w) for _, _, w in edges)/(4*arb.pi()*factorial(kernel["ORDER"]+1))

    def weights(k, w, d):
        return [k, w, (1-t)*d, t*d, -d]

    def tail(bg, coeff):
        size = arb(0)
        for (g, a, _), charges in zip(bg, fields, strict=True):
            size += abs(2*coeff[0]*a) + abs(coeff[1]*g)
            size += sum(abs(c*q) for c, q in zip(coeff[2:], charges, strict=True))
        return prefactor*size

    bounds = (".01180", ".03252", ".02986")
    for (i, j), direction, bound in zip(combinations(range(3), 2), DIRECTIONS, bounds, strict=True):
        p, k, u, h, v = map(arb, direction)
        left = weights(k+h, v, p+u)
        right = weights(-h, -v, -u)
        coeff = arb_mat([left])*series[i] + arb_mat([right])*series[j]
        error = tail(rows[i], left) + tail(rows[j], right)
        lower = (p*intercept - k + (coeff*moments)[0, 0] - error).lower()
        # 精确匹配引力时，只余 (p+u)D_i-uD_j；两处可独自选能标。
        gap = lower/((abs(p+u)+abs(u))*gauge_norm)
        print(HEAVY[i], HEAVY[j], "共同分布的规范差下界", gap.str(12), flush=True)
        assert gap > arb(bound)

    # 第一和第三个重谱单靠规范差就能区分，无需用引力条件。
    pure = weights(0, 0, 1)
    diff = arb_mat([pure])*(series[0]-series[2])
    upper = ((diff*moments)[0, 0] + tail(rows[0], pure) + tail(rows[2], pure)).upper()
    print("不约束引力时第一与第三重谱的 D 差上界", upper.str(12))
    assert upper < arb("-.02435")

    solve_witness(load, shell, 323)


def solve_witness(load, shell, line, *, layout=None, energy=10**12, minimum="11.58"):
    if layout is None:
        layout = json.loads((ROOT / f"ch39-fifth-shell-common-{line}.json").read_text())
    labels = [arb(k)/11000 for row in layout["base"] for k in row]
    bases = []
    for k in layout[f"line_{line}"]:
        fixed = labels.copy()
        fixed[line] = arb(k)/11000
        bases.append(fixed)
    extra = [(who, i, arb(k)/11000) for who, i, k in layout["extra"]]
    load("one-background-specific-line")["solve_layout"](
        bases, extra, line, [minimum]*3, [shell]*3, HEAVY, False, reference_gev=energy)


if __name__ == "__main__":
    verify()
