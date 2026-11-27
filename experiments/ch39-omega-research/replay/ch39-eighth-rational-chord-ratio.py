"""在同一有理球面位置按弦距比和按半角坐标分别读取第八壳第二枢轴。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb


CENTERS = (
    "55.098018319140046",
    "76.67695591387131",
    "67.65284285664114",
    "76.05978748157519",
    "0.6818144605343499",
    "75.99829769093567",
    "63.04633487955998",
    "1.9909487309502854",
    "9.555005605006956",
    "0.4196856574038982",
    "71.4589352655587",
    "1.7447316732959",
    "7.472079483277022",
    "64.25583001419216",
    "76.70899196638433",
    "0.5092289398242106",
    "52.98628865537255",
    "52.355117811220204",
    "27.382879482573685",
)
ACTIVE = (15, 16, 14, 0, 8)


def squared_distance(a, b):
    return sum((x - y)**2 for x, y in zip(a, b, strict=True))


def verify(root):
    a, b = (23, 13, 13), (25, 11, 11)
    p = tuple(Fraction(k, 171) for k in (4195, 1969, 1969))
    assert sum(x*x for x in p) == 867
    assert squared_distance(p, a) / squared_distance(p, b) == Fraction(289, 32)
    r = Fraction(3, 4)
    t = 17*r / (24 - 7*r)
    u, v = (1 - t*t)/(1 + t*t), 2*t/(1 + t*t)
    assert (t, u, v) == (Fraction(17,25), Fraction(168,457), Fraction(425,457))
    assert (1-u)/(1-v) == Fraction(289, 32)
    another_u, another_v = (1-r*r)/(1+r*r), 2*r/(1+r*r)
    assert (another_u, another_v) == (Fraction(7,25), Fraction(24,25))
    print("有理位置与两端平方弦距的比", p, Fraction(289,32), flush=True)

    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    arc = lib("eighth-shortest-arc-pivot-support")
    support = lib("cross-shell-pivot-support")
    tables = support["projected_tables"](lib("fourth-third-pivot-global-kinetic-obstruction"),
                                         root, arc["NUMBER"], None)
    rows = arc["polynomials"](root, include_single_edge=False)[0][1]
    rational = lib("eighth-rational-arc-signs")
    for parameter, sign in ((r, 1), (t, -1)):
        weights = rational["rational_weights"](rows, rational["interval"](parameter))
        lo, hi = map(int, arc["score"](tables, weights))
        assert (lo >= 2_834_300 if sign == 1 else hi < 0)
        print("半角参数与指标的有向包围", parameter, lo, hi, arc["SCALE"], flush=True)
    lib("eighth-arc-positive-probe")["verify"](
        root, parameter_case=(arb(17) / 25, CENTERS, ACTIVE)
    )


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
