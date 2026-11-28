"""固定第五壳及参考分离方向，以同一个壳矩区间检验六个留出重谱。"""
import json
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
REFERENCE = 2_345_674
HOLDOUT = tuple(range(2_345_676, 2_345_687, 2))


def shell_moments(edges, kernel, scale):
    """每个单调基函数在每条线上取全区间，再以有符号权求和。"""
    moments = [arb(0) for _ in range(kernel["ORDER"] + 4)]
    for left, right, weight in edges:
        assert 0 < left and (left < right or left == right) and right <= 1
        lo = kernel["basis"](scale * left).entries()
        hi = kernel["basis"](scale * right).entries()
        for j, (a, b) in enumerate(zip(lo, hi, strict=True)):
            moments[j] += weight * a.union(b)
    return arb_mat([[entry] for entry in moments])


def verify():
    exchange = run_path(str(ROOT / "ch39-fourth-fifth-shell-exchange.py"))
    four = run_path(str(ROOT / "ch39-four-background-crossing.py"))
    load = four["load"]
    core = load("shared-large-scale-gap")
    context = core["load"]()
    fields, kernel = context[1], context[-2]
    edges = load("shared-twelve-node-gauge-gap")["paired_lines"](
        ROOT, load("native-rh-ray-pvm")["shell_weights"], exchange["FIFTH_SHELL"])
    assert len(edges) == 361
    scale = arb(7) / 1000
    floor = scale * min(a.lower() for a, _, _ in edges)
    moments = shell_moments(edges, kernel, scale)
    direction = exchange["SEPARATOR"]
    ky, wy, u, v = [arb(n) / 10**11 for n in direction]
    projection = arb_mat([[ky, wy, u, v, -u-v]])
    numbers = (REFERENCE, *HOLDOUT)
    _, _, spectra = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT, numbers, tight_neutral=True)
    lower_bounds = dict(zip(HOLDOUT, [
        ("-.3771", "-.1895"), (".1678", ".2091"), ("-.4015", ".9926"),
        (".1356", ".2423"), ("-.4252", "-.2237"), (".1626", ".1821"),
    ], strict=True))
    transferred, excluded = [], []
    anchor = None
    for number, rows in zip(numbers, spectra, strict=True):
        assert max(z.upper() for _, _, z in rows) < min(arb(1), floor)
        series = kernel["series"](rows, fields)
        single = ([(edges, series, rows)], *context[1:])
        groups, constant, norm = core["project_direction"](single, direction, kernel_ratio=1)
        direct = four["lower_projection"](single, groups, constant, 0)
        coefficients = projection * series
        if number == REFERENCE:
            assert direct / norm > arb(".2359")
            anchor = (direct, constant, coefficients)
            continue
        reference, old_constant, old_coefficients = anchor
        # constant 已扣本端的尾界。先加回参考端扣除的常数，再换到本端。
        delta = ((coefficients - old_coefficients) * moments)[0, 0]
        transfer = reference + (constant - old_constant) + delta.lower()
        expected_transfer, expected_direct = lower_bounds[number]
        assert transfer/norm > arb(expected_transfer) and direct/norm > arb(expected_direct)
        if transfer > 0:
            transferred.append(number)
        if direct > 0:
            excluded.append(number)
        print(number, "壳矩迁移", (transfer/norm).str(12),
              "逐线", (direct/norm).str(12), flush=True)
    assert transferred == [2_345_678, 2_345_682, 2_345_686]
    assert excluded == [2_345_678, 2_345_680, 2_345_682, 2_345_686]

    layouts = json.loads((ROOT / "ch39-fifth-shell-holdout-layouts.json").read_text())
    assert [entry["heavy"] for entry in layouts] == [n for n in HOLDOUT if n not in excluded]
    solve = load("one-background-specific-line")["solve_layout"]
    for entry in layouts:
        labels = [arb(k)/11000 for row in entry["base"] for k in row]
        extra = [(None, i, arb(k)/11000) for i, k in entry["extra"]]
        print("单独构造的留出重谱", entry["heavy"], flush=True)
        solve([labels], extra, 0, ["11.6"], [exchange["FIFTH_SHELL"]],
              [entry["heavy"]], False, reference_gev=10**12)


if __name__ == "__main__":
    verify()
