"""第八壳第二枢轴的 A_c≥2 支撑判据在两种弧读出之间的临界。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    support = lib("cross-shell-pivot-support")
    arc = lib("eighth-shortest-arc-pivot-support")
    rational = lib("eighth-rational-arc-signs")
    rows = arc["polynomials"](root, include_single_edge=False)[0][1]
    raw = support["projected_tables"](core, root, arc["NUMBER"], None)
    tables = []
    for (offset, bands, scale), d in zip(raw, support["DIRECTIONS"], strict=True):
        assert d[3] <= 0
        tables.append((offset - d[3]*arb(2)*core["SCALE"]/10, bands, scale))
    centers = (
        "75.17165038503167", "5.314191467656551", "77.09181518401721",
        "64.41308469069699", "7.095530004705824", "14.11513983226248",
        "0.26700844227013637", "39.00753982114472", "7.104499850307511",
        "0.27212836556978104", "77.90383494222813", "36.009425529046055",
        "0.33545827803194467", "41.92832563195316", "11.225309946451707",
        "0.2703893698985129", "55.92250846917611", "34.785371043130176",
        "6.029844966359748",
    )
    parameter = arb(".734") + arb(0, arb("5e-11").upper())
    lib("eighth-arc-positive-probe")["verify"](root, parameter_case=(parameter, centers, (15,10,16,8,7)))
    for a, b, sign in (("0.668", "0.73615841", -1), ("0.73615844", "0.75", 1)):
        rational["check_sign"](arc, tables, rows, a, b, sign)
    rational["crossing"](arc, tables, rows, "0.73615841", "0.73615844", (0,7,16), 1)
    for t in ("0.73615841","0.73615844"):
        weights = rational["rational_weights"](rows, arb(t))
        bounds = arc["score"](tables, weights)
        print("临界窗端点、平移指标的整数包围（单位一千万分之一）", t, *map(int,bounds), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
