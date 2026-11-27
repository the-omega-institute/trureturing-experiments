"""在十三个实际壳上直接读末端口八站的荷五传播和等模相干端点。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, acb, ctx


def verify(root):
    graph = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    sites = run_path(str(root / "ch39-terminal-port-support.py"))["verify"](root)
    pairs = list(combinations(sites,2))
    expected = [(3,6),(3,5),(9,13),(9,13),(9,13),(9,13),(5,6),
                (5,6),(5,6),(3,6),(3,6),(3,5),(3,5)]
    gaps = []
    for shell,pair in enumerate(expected, start=7):
        ctx.prec = 256
        height = acb.zeta_zero(2_306_888_100_085+shell).imag
        green = graph["laplacian"](5,height).inv()
        values = {edge:abs(green[edge]) for edge in pairs}
        value = values[pair]
        assert all(x.is_finite() for x in values.values())
        assert all(value < x for edge,x in values.items() if edge != pair)
        # 固定一个端点在 6，另一个端点以等模叠加在 9、13，对相位取最小。
        difference = (abs(green[6,9])-abs(green[6,13]))/2
        assert difference.is_finite()
        assert difference>0 if shell<=11 or shell>=18 else difference<0
        gaps.append(abs(difference))
        print("壳号、最小交换站对、点流振幅模、等模端点的最小振幅",
              shell,pair,value.str(13),gaps[-1].str(13))
    assert arb("4.36e-7") < gaps[4] < arb("4.37e-7")
    assert all(gaps[4] < g for j,g in enumerate(gaps) if j!=4)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
