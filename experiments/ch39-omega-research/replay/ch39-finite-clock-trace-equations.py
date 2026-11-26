"""从十九点壳的闭游走列出单背景谱反演的有限 Laurent 方程与代数次数上界。"""
from collections import defaultdict
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def clock_equations(shell, delta, charge, longest):
    neighbors = [[] for _ in range(19)]
    diagonal = [arb(0) for _ in neighbors]
    pole, tuned = shell["soft"]["pole"], shell["soft"]["star"]["orbit"]
    for i, j, distance, unit in shell["geo"]["edges"]:
        spoke = pole in (i, j)
        weight = 1 - delta if spoke and (i in tuned or j in tuned) else arb(1)
        diagonal[i] += weight
        diagonal[j] += weight
        value = -acb(weight) if spoke else -weight * unit**charge
        a, b = charge * int(distance == 724), charge * int(distance == 804)
        assert not spoke or a == b == 0
        neighbors[i].append((j, a, b, value))
        neighbors[j].append((i, -a, -b, value.conjugate()))
    for i, d in enumerate(diagonal):
        neighbors[i].append((i, 0, 0, d))
    traces = [defaultdict(acb) for _ in range(longest + 1)]
    for start in range(19):
        paths = {(start, 0, 0): acb(1)}
        for length in range(1, longest + 1):
            steps = defaultdict(acb)
            for (i, a, b), coefficient in paths.items():
                for j, da, db, value in neighbors[i]:
                    steps[j, a + da, b + db] += coefficient * value
            for (i, a, b), coefficient in steps.items():
                if i == start:
                    traces[length][a, b] += coefficient
            paths = steps
    return traces


def hexagon(m):
    return {(a, b) for a in range(-m, m + 1) for b in range(-m, m + 1) if abs(a+b) <= m}


def verify():
    shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
    ctx.prec = 256
    h = acb.zeta_zero(2_306_888_100_085).imag
    z, w = (shell["geo"]["phases"](h)[d] for d in (724, 804))
    delta = arb(23)/25 * arb(2)**-20 / (h*h + shell["soft"]["point"]["study"]["offset"])
    for charge, last in ((1, 4), (3, 3)):
        traces = clock_equations(shell, delta, charge, last)
        native = shell["laplacian"](charge, h)
        for order in range(1, last + 1):
            radius = 0 if order < 3 else order - 2
            expected = {(charge*a, charge*b) for a, b in hexagon(radius)}
            assert set(traces[order]) == expected
            value = sum(c * z**a * w**b for (a, b), c in traces[order].items())
            matrix = native**order
            assert (value - sum(matrix[j, j] for j in range(19))).contains(0)
            print("荷、阶数、闭游走的不同指数数目", charge, order, len(expected))
    # mH 的面积为 3m²，两方程的一二倍与一三倍六边形分别给 12 与 18。
    for m, target in ((2, 12), (3, 18)):
        area_difference = 3 * ((m+1)**2 - 1 - m*m)
        assert area_difference == target
        print("第二个六边形的倍率、孤立单根个数上界", m, target)


if __name__ == "__main__":
    verify()
