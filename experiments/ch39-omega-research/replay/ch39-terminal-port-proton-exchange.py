"""末端口八站上荷五原生传播的正电流锐界与同模相干抵消。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    selected = run_path(str(root / "ch39-terminal-port-support.py"))["verify"](root)
    exchange = run_path(str(root / "ch39-shell-current-word-floor.py"))
    geometry = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))["geo"]
    points = exchange["points"]
    green = exchange["green"]
    assert list(points) == list(geometry["sites"])
    assert selected == (0, 2, 3, 5, 6, 9, 11, 13)
    assert all(green[i, j].is_finite() for i in selected for j in selected)
    assert arb(".02424488184856") < green[9, 13].real < arb(".02424488184857")
    magnitude = abs(green[9, 13])
    assert arb(".02445238661631") < magnitude < arb(".02445238661632")
    for first, second in ((9, 13), (13, 9)):
        # 把源概率单纯形切为 a_first >= a_second 的两部分。新产生的唯一顶点为中点。
        direction = green[first, second].conjugate() / abs(green[first, second])
        margins = []
        for i in selected:
            if i == second:
                continue
            for j in selected:
                if (i, j) == (first, second):
                    continue
                margins.append((direction * green[i, j]).real.lower())
        for j in selected:
            midpoint = (green[first, j] + green[second, j]) / 2
            margins.append((direction * midpoint).real.lower())
        assert all(value.is_finite() and value > arb(".029") for value in margins)
        print("取等点、其余顶点最小有向投影", first, second, min(margins).str(15))
    print("正端点电流的最小传播模", magnitude.str(16))

    # 反射固定极点，保留点积而翻转有向三重积。在已固定的极点规范中，
    # 相应反酉对称逐边给出 P L5 P^T = conj(L5)，取逆便得到 G5 的恒等式。
    pole = exchange["pole"]
    reflection = [points.index((y, x, z)) for x, y, z in points]
    assert {reflection[i] for i in selected} == set(selected)
    assert reflection[6] == 6 and reflection[9] == 13 and reflection[13] == 9
    dot, det = geometry["dot"], geometry["det"]
    for i, j in combinations(range(len(points)), 2):
        x, y = points[i], points[j]
        rx, ry = points[reflection[i]], points[reflection[j]]
        assert dot(pole, x) == dot(pole, rx) and dot(y, pole) == dot(ry, pole)
        assert dot(x, y) == dot(rx, ry) and det(pole, x, y) == -det(pole, rx, ry)
    link = green[6, 9]
    assert abs(link) > arb(".029")
    phase = link / green[6, 13]
    # 模一来自上述精确反射恒等式。这里的区间重叠只核对实现。
    assert abs(phase).contains(1)
    assert ((link - phase * green[6, 13]) / 2).contains(0)
    print("同模端点的非零分支", abs(link).str(16))
    return magnitude


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
