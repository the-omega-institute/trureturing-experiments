"""第八、九壳第二枢轴在 t=3/4 的谱权与单层响应交叉检验。"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    arc = lib("eighth-shortest-arc-pivot-support")
    rational = lib("eighth-rational-arc-signs")
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    support = lib("cross-shell-pivot-support")
    weights, tables = {}, {}
    for shell in (8, 9):
        n = 2_306_888_100_085 + shell
        rows, _ = arc["polynomials"](root, number=n, include_single_edge=False)
        w = rational["rational_weights"](rows[1], arb(3)/4)
        assert len(w) == 19 and sum(w).contains(1)
        order = sorted(range(19), key=lambda i: abs(w[i]).mid(), reverse=True)
        assert order[:3] == [0, 16, 7]
        assert min(abs(w[i]).lower() for i in order[:3]) > max(abs(w[i]).upper() for i in order[3:])
        tables[shell] = support["projected_tables"](core, root, n, None)
        weights[shell] = w

    def read(w, table):
        box = arc["score"](table, w)
        assert box is not None
        return tuple(Fraction(int(x), arc["SCALE"]) for x in box)

    for shell, w in weights.items():
        for response, table in tables.items():
            lo, hi = read(w, table)
            print("谱权壳、响应壳、S3 有理闭区间", shell, response, lo, hi, flush=True)
            assert (lo > 0 if shell == 8 else hi < 0)

    old, new = weights[8], weights[9]
    delta = [b-a for a, b in zip(old, new, strict=True)]
    assert all(x > 0 or x < 0 for x in delta)
    flips = []
    tested = 0
    # 两项按实际差额相向移动，取一个层到达新值时的最大步；权和仍为一。
    for i, j in combinations(range(19), 2):
        if delta[i]*delta[j] > 0:
            continue
        up, down = (i, j) if delta[i] > 0 else (j, i)
        amount = delta[up].min(-delta[down])
        w = list(old)
        w[up] += amount
        w[down] -= amount
        assert sum(w).contains(1)
        left, right = [read(w, tables[k]) for k in (8, 9)]
        tested += 1
        if left[1] < 0 and right[1] < 0:
            flips.append((up+1, down+1))
            print("增权层、减权层、转移动量、两表的 S3 区间",
                  up+1, down+1, amount.str(10), left, right, flush=True)
        else:
            assert left[0] > 0 and right[0] > 0
    assert set(flips) == {(10, 17), (16, 17)}
    print("受检同向补偿层对数、其中双表变号数", tested, len(flips), flush=True)

    # 三个最大绝对权之外的同号系数仅以总量进入检验，两个受益层在相同步长时等效。
    for layer in (10, 16):
        w = list(old)
        assert delta[layer-1] > arb(3)/4 and -delta[16] > arb(3)/4
        w[layer-1] += arb(3)/4
        w[16] -= arb(3)/4
        assert min(abs(w[i]).lower() for i in (0,7,16)) > max(
            abs(x).upper() for i,x in enumerate(w) if i not in (0,7,16)
        )
        results = [read(w, tables[k]) for k in (8,9)]
        assert all(hi < Fraction(-1779,10000) for lo, hi in results)
        print("只转移 3/4、增权层、两表区间", layer, *results, flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
