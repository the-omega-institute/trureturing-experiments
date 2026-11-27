"""第九壳第二枢轴在三个预登记有理读出上的首次留出检验。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path


NUMBER = 2_306_888_100_094
PREDICTIONS = ((3, 4, 5, 1), (5, 12, 13, -1), (7, 24, 25, 1))


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    arc = lib("eighth-shortest-arc-pivot-support")
    rational = lib("eighth-rational-arc-signs")
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    tables = lib("cross-shell-pivot-support")["projected_tables"](core, root, NUMBER, None)
    rows, _ = arc["polynomials"](root, number=NUMBER, include_single_edge=False)
    scale = arc["SCALE"]
    for a, b, c, expected in PREDICTIONS:
        assert a*a + b*b == c*c
        t = Fraction(b, c+a)
        weights = rational["rational_weights"](rows[1], rational["interval"](t))
        assert sum(weights).contains(1)
        box = arc["score"](tables, weights)
        if box is None:
            print(t, "前三层未分离，未定", flush=True)
            continue
        lo, hi = map(int, box)
        actual = 1 if lo > 0 else -1 if hi < 0 else 0
        outcome = "未定" if actual == 0 else "命中" if actual == expected else "失配"
        print("半角参数、支撑指标有理闭区间、预登记符号、结果",
              t, Fraction(lo, scale), Fraction(hi, scale), expected, outcome, flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
