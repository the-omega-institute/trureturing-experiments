"""联合编号起的八项规范不变谱读数和有界频率反演的区间前提。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    native = run_path(str(root/"ch39-unlabeled-native-connection.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    faces = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 320
    sites = shell["geo"]["sites"]
    distances = [[sum((a-b)**2 for a, b in zip(x, y, strict=True))
                  for x in sites] for y in sites]
    pos, mirrored = (native["reconstruct"](distances, i) for i in (0, 1))
    # 保持空间取向和两个钟，只把软化的远点轨道换到其镜像。
    mirror_edges = {pair: (pos[-1][pair][0], clock, direction, soft)
                    for pair, (_, clock, direction, soft) in mirrored[-1].items()}
    neg = (*mirrored[:-1], mirror_edges)
    start = 2_306_888_100_085
    heights = [z.imag for z in acb.zeta_zeros(start, 8)]
    assert all(a < b for a, b in zip(heights, heights[1:]))
    width = heights[-1]-heights[0]
    assert arb("1.6073") < width < arb("1.6074") < 2
    x, y = (faces["unit"](v) for v in ((510, 1), (1173, -2)))
    assert not x.imag.contains(0) and not y.imag.contains(0)
    assert 0 < shell["geo"]["a"] < shell["geo"]["b"] < arb.pi()/2
    readings = []
    for h in heights:
        phases = native["clocks"](pos, h)
        delta = arb(23)/25 * arb(2)**-20 / (h*h+shell["soft"]["point"]["study"]["offset"])
        difference = native["matrix"](pos, 1, phases, delta)**3 - native["matrix"](neg, 1, phases, delta)**3
        actual = sum(difference[i, i] for i in range(19)) / delta
        analytic = -36*(x.imag*phases[1].imag + y.imag*phases[2].imag)
        assert actual.imag.contains(0) and (actual.real-analytic).contains(0)
        assert actual.real.rad() < arb("1e-45")
        readings.append(actual.real.str(14))
    print("首末实际零点的高度差、八个荷一三阶归一化谱读数", width.str(17), readings)


if __name__ == "__main__":
    verify()
