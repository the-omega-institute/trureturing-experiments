"""第八壳的第一枢轴端点排除与第二枢轴五式精确根。"""
from argparse import ArgumentParser
from pathlib import Path
from runpy import run_path

from flint import arb


NUMBER = 2_306_888_100_093
CENTERS = (
    "2.4167061617364407e-34", "2.124407583177014e-31", "1.5660541858747768e-31",
    "2.3594085647379127e-31", "4.555985718151624e-32", "1.33849992840502e-31",
    "3.15787995747478e-12", "3.694474989801832e-32", "4.332646366234258e-31",
    "0.31017631035415977", "0.001258452850144041", "2.4327270981468317e-34",
    "4.435664259549921e-31", "3.0697820753938366e-31", "2.1750770588544744e-18",
    "3.399594768640684e-11", "4.51898583756821e-18", "0.08914044020202493",
    "0.08649908740640672",
)


def first(root):
    exclude = run_path(str(root / "ch39-fifth-arbitrary-pivot-obstructions.py"))["verify"]
    return exclude(root, number=NUMBER, pivots_to_check=(1,), dominant=16)


def second(root):
    solve = run_path(str(root / "ch39-fifth-third-pivot-em-colour-path.py"))["solve"]
    active = (17, 10, 1, 18, 16, 11)
    data, rows, weights, sites, boxes, _, lib = solve(
        root, CENTERS, tuple(i - 1 for i in active), number=NUMBER, pivot=2
    )
    gauge = lib("all-field-signed-crossing")
    direct = lib("full-neutral-running")["direct"]
    rule = gauge["angular_rule"](16, 45)
    norm = sum(abs(arb(coeff[2])) for _, _, coeff in rows) / (4 * arb.pi())
    motion = sum(
        abs(w) * abs((box / x).log()) for w, box, x in zip(weights, boxes, sites, strict=True)
    )
    err = arb(0, (norm * motion).upper())
    colour = []
    for energy in ("1.176397e14", "560000"):
        value = direct(data, sites, weights, rows, energy, gauge, rule)[2] + err
        print("能标 GeV、第二枢轴精确根的色逆耦合", energy, value.str(13), flush=True)
        colour.append(value)
    assert arb("35.783") < colour[0] < arb("35.785")
    assert arb("-95.274") < colour[1] < arb("-95.272")
    return colour


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--pivot", type=int, choices=(1, 2))
    chosen = parser.parse_args().pivot
    root = Path(__file__).resolve().parent
    if chosen in (None, 1):
        first(root)
    if chosen in (None, 2):
        second(root)
