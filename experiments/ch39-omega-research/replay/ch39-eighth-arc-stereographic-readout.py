"""第八壳同一球面点上，两种同步规则的色逆耦合界。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    arc = lib("eighth-shortest-arc-pivot-support")
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    tables = lib("cross-shell-pivot-support")["projected_tables"](core, root, arc["NUMBER"], None)
    rows = arc["polynomials"](root, include_single_edge=False)[0][1]
    rational_weights = lib("eighth-rational-arc-signs")["rational_weights"]

    # 球面弦长平方为 12，半径平方为 867，角 alpha 为三条最短弧的圆心角。
    cos_alpha = 1 - arb(12) / (2 * 867)
    tan_half_squared = (1 - cos_alpha) / (1 + cos_alpha)
    assert (tan_half_squared * 288).contains(1)
    cos_half = ((1 + cos_alpha) / 2).sqrt()
    cos_quarter = ((1 + cos_half) / 2).sqrt()
    cos_point = 4 * cos_quarter**3 - 3 * cos_quarter
    t_projection = ((1 - cos_point) / ((1 + cos_point) * tan_half_squared)).sqrt()

    # 原角度同步在量子态圆上到达 3pi/8，用两次半角公式给出 t=tan(3pi/16)。
    cos_3pi8 = ((1 - 1 / arb(2).sqrt()) / 2).sqrt()
    t_angle = ((1 - cos_3pi8) / (1 + cos_3pi8)).sqrt()
    assert arb(".66") < t_angle < arb(".68")
    assert arb(".749") < t_projection < arb(".751")

    for label, parameter, positive in (
        ("角度分数同步", t_angle, False),
        ("归一化半角坐标同步", t_projection, True),
    ):
        weights = rational_weights(rows, parameter)
        assert len(weights) == 19 and sum(weights).contains(1)
        lo, hi = arc["score"](tables, weights)
        if positive:
            assert lo >= 2_825_459
        else:
            assert hi <= -8_300_460
        print(label, "参数", parameter, "S3 有向包围", int(lo)/arc["SCALE"],
              int(hi)/arc["SCALE"], flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
